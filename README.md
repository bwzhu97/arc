# Numerical experiments for Riemannian ARC subproblem solvers

This repository contains the Matlab codes to reproduce our numerical experiments comparing different subproblem solver implementations of Riemannian Adaptive Regularization with Cubics (Riemannian ARC).

The main script is compare_solvers_on_problems.m.

It is necessary to download Manopt (and to add it to Matlab's path) to run these experiments. Preferably obtain the latest version from https://github.com/NicolasBoumal/manopt. Please run with a version of Manopt dating from Oct. 4, 2018 or more recent.

The code pdf_print_code.m produces PDFs from Matlab figures. It relies on the utility pdfcrop being installed. If this is an issue, simply remove the call to pdf_print_code in the main script.

Authors:
Bryan Zhu and Nicolas Boumal (based on code authored by Naman Agarwal, Nicolas Boumal, Brian Bullins, Coralia Cartis at https://github.com/NicolasBoumal/arc)

May 2, 2019
