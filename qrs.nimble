# Package

version       = "0.1.0"
author        = "Gennaro Calendo"
description   = "Quick coverage statistics from per-base BED files"
license       = "MIT"
srcDir        = "src"
bin           = @["qrs"]


# Dependencies

requires "nim >= 2.2.10"
requires "lapper >= 0.1.8"
requires "hts >= 0.3.31"
requires "cligen >= 1.9.6"