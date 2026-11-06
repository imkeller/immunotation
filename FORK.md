# About this fork

This is [Bioconductor's immunotation](https://bioconductor.org/packages/immunotation) with fixes for the
Allele Frequency Net Database (AFND) moving its pages. Install it with:

```r
pak::pkg_install("fmobegi/immunotation")
```

## What differs from Bioconductor

- Every AFND URL points to `https://www.allelefrequencies.net/afnd-legacy/`. The old addresses
  (`hla6006a.asp`, `hla6006c.asp`, `hla6003x.asp`, `pop6001c.asp`) return 404, which made the package
  impossible to install from source (always the case on Linux).
- `valid_geographics`, which is downloaded while the package loads, falls back to empty lists if AFND is
  unreachable, so the package still installs.
- Version: Bioconductor's version plus a fourth number (e.g. `1.21.0.1`), so it is always newer than the
  Bioconductor copy it is based on.

## Pulling in Bioconductor updates

The `bioc` remote is Bioconductor's own git server (the GitHub repository `imkeller/immunotation` stopped
at 1.7.0). In a fresh clone, add it once:

```bash
git remote add bioc https://git.bioconductor.org/packages/immunotation
```

Then, whenever Bioconductor releases a new version:

```bash
git fetch bioc
git merge bioc/devel
```

If the merge conflicts, keep the `afnd-legacy` URLs and the `tryCatch()` around `get_valid_geographics()`,
and set `Version:` to Bioconductor's new version plus `.1`. Check with `git diff bioc/devel`: only those
fixes should remain.
