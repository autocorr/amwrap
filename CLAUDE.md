# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Install (also compiles AM from C source):**
```bash
pip install .
```

**Run all tests:**
```bash
python -m pytest amwrap/tests/test_core.py
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

The entire Python module lives in `amwrap/__init__.py`. There is no subpackage structure. `amwrap/plotting.py` is a standalone utility for figure output and is not part of the core model pipeline.

### Module-level initialization

At import time, `amwrap/__init__.py` does several things that affect all subsequent calls:

- **`AmExecutable`** probes the system PATH for `am` and `am-serial` binaries; falls back to the bundled binaries in `amwrap/bin/`. Two singletons are created: `AM_PARALLEL` (OpenMP build) and `AM_SERIAL` (single-threaded build).
- **`ENV`** is a frozen copy of `os.environ` with `AM_CACHE_PATH` set to `/dev/shm/` when that directory exists. This dict is passed to all subprocess calls. It must not be mutated; use `{**ENV, ...}` to override keys per-call.
- **`CLIMATOLOGIES`** eagerly loads all six standard atmospheric profiles from `amwrap/climatology/*.dat` into `Climatology` instances.

### Model pipeline

`Model` holds atmospheric configuration as `astropy.units.Quantity` arrays (pressure, temperature, per-species mixing ratios, optional cloud layers). The `.config_text` property renders this into an AM configuration file format, which is passed to AM via stdin. AM's output is the NumPy binary (`.npy`) format streamed to stdout; `_parse_output()` loads it with `np.load(BytesIO(...))` and packages it into a `pandas.DataFrame`.

`Model.run(parallel=False, cache_dir=None)`:
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
