# NICE

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://quantumelephant.github.io/NICE.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://quantumelephant.github.io/NICE.jl/dev/)
[![Build Status](https://github.com/quantumelephant/NICE.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/quantumelephant/NICE.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/quantumelephant/NICE.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/quantumelephant/NICE.jl)
# NICE
NICE is a package for determining solving simultaneous equilibria problems with networks of complex chemical reactions. The name comes because the main approach is **N**et-event kinetic Monte Carlo, applied to fill in an **ICE** ([Initial, Charge, Equilibrium](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Equilibria/Le_Chateliers_Principle/Ice_Tables) table.

The scope of `NICE` also supports direct solution of the nonlinear equations associated with simultaneous equilibria and (normal) kinetic Monte Carlo. While some innovations are present in the net-event kinetic Monte Carlo code, the intention of NICE is to support studies of relatively complicated but small systems: sophisticated parallelization techniques that would be needed to treat very large-scale problems are not included.
