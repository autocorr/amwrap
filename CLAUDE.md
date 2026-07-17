# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Install (also compiles AM from C source):**
```bash
pip install .
```

**Run all tests:**
```bash
python -m pytest amwrap/tests/
```

**Run a single test:**
```bash
python -m pytest amwrap/tests/test_core.py::TestModel::test_run
```

**Rebuild AM binaries only** (after modifying C source in `src/am-14.0/`):
```bash
cd src/am-14.0
make serial && mv am ../../amwrap/bin/am-serial && make clean
make am    && mv am ../../amwrap/bin/am         && make clean
cd ../..
```

**Clean build artifacts:**
```bash
make clean
```

## Architecture

AM is invoked as a subprocess, **by design**. Direct C bindings (shared library + ctypes) were evaluated and rejected: AM keeps mutable global state (the `output[]` table, errlog, kcache), reads `AM_CACHE_PATH` once per process into a static (which would break the per-call `cache_dir` override), a C crash would take down the interpreter, and the vendored source currently upgrades patch-free. Process isolation is also what makes per-worker `multiprocessing` parallelism safe. Do not propose in-process bindings without revisiting these constraints.

The package is a facade over five submodules; `amwrap/__init__.py` re-exports the public API listed in `__all__`:

- `amwrap/constants.py` — physical constants and `MOD_DIR`.
- `amwrap/thermo.py` — free thermodynamic/interpolation helpers.
- `amwrap/climatology.py` — `Climatology` and the `CLIMATOLOGIES` dict.
- `amwrap/driver.py` — everything that knows AM is a subprocess: `AmExecutable` discovery, the frozen `ENV`, `get_executable()`, and `run_am()` (the injectable execution seam).
- `amwrap/model.py` — `Model`: validation, `config_text` rendering, `_parse_output()`.

`amwrap/plotting.py` is a standalone utility for figure output and is not part of the core model pipeline. Import discipline: never `from .driver import AM_PARALLEL` or `from .climatology import CLIMATOLOGIES` — always access them as module attributes, or the lazy loading below is defeated.

### Module-level initialization

- **`AM_PARALLEL`** (OpenMP build) and **`AM_SERIAL`** (single-threaded build) singletons probe the system PATH for `am` and `am-serial`, falling back to the bundled binaries in `amwrap/bin/`. Probing is **lazy** (PEP 562 module `__getattr__`, cached into module globals): it runs on first access of any `AM_*` name or first `Model.run()`, not at import. The missing-binary and version-mismatch `UserWarning`s fire at that point.
- **`CLIMATOLOGIES`** likewise lazily loads all six standard atmospheric profiles from `amwrap/climatology/*.dat` on first access.
- **`ENV`** is built eagerly at import: a frozen copy of `os.environ` with `AM_CACHE_PATH` set to `/dev/shm/` when that directory exists. This dict is passed to all subprocess calls. It must not be mutated; use `{**ENV, ...}` to override keys per-call.

### Model pipeline

`Model` holds atmospheric configuration as `astropy.units.Quantity` arrays (pressure, temperature, per-species mixing ratios, optional cloud layers). The `.config_text` property renders this into an AM configuration file format, which is passed to AM via stdin. AM's output is the NumPy binary (`.npy`) format streamed to stdout; `_parse_output()` loads it with `np.load(BytesIO(...))` and packages it into a `pandas.DataFrame`.

`Model.run(parallel=False, cache_dir=None)` delegates execution to `driver.run_am()`, which takes the executable and env as parameters (inject a stub there to test the run path without a compiled AM; see `tests/test_driver.py`):
- `parallel=True` selects the OpenMP build (`am`); `False` selects the serial build (`am-serial`).
- `cache_dir` overrides `AM_CACHE_PATH` for that subprocess call. Use this when running many models in parallel via `multiprocessing` — AM's per-bucket file locking is not exclusive on POSIX (`rename()` replaces atomically), which causes `rename_with_retry()` to sleep up to 39 s per collision under contention. Giving each worker process its own cache subdirectory eliminates this:

```python
def worker(model_kwargs):
    cache = Path(f"/dev/shm/amcache_{os.getpid()}")
    cache.mkdir(exist_ok=True)
    return Model(**model_kwargs).run(cache_dir=cache)
```

### AM disk cache (`dcache.c`)

Controlled by two environment variables:
- `AM_CACHE_PATH`: directory for cache files (default: none/disabled).
- `AM_CACHE_HASH_MODULUS`: number of hash buckets (default: 1021). Cache files are named `am_<hashval>_<slot>` with up to 4 slots per bucket.

The cache keys on `(df, P, vmr, T, lineshape, k_type)` with exact floating-point comparison. Models with identical parameters share cache entries; models varying only in cloud liquid may still share all non-cloud absorption coefficients.

### Building AM

`setup.py` drives the build at install time: it `cd`s into `src/am-14.0/`, runs `make serial` then `make am`, and moves the resulting binaries to `amwrap/bin/`. The two targets differ only in whether OpenMP is enabled. The C source is AM 14.0 by Scott Paine (SAO); do not modify it without understanding the full spectroscopic model.
