#!/bin/bash
Rscript ~/RDM/CARDAMOM/example_postprocessing/Create_latex_C_budget_script_v2_TG.r
pdflatex ./latex_C_budget_TGtest.tex 
pdftoppm -png latex_C_budget_TGtest.pdf latex_C_budget_TGtest
rm latex_C_budget_TGtest.pdf