# PETRELS-ADMM: Robust Subspace Tracking with Missing Data and Outliers

We propose a novel algorithm called PETRELS-ADMM to deal with subspace tracking in the presence of outliers and missing data. The proposed approach consists
of two main stages: outlier rejection and subspace estimation. Particularly, we first use ADMM solver for detecting outliers corrupted in data in an efficient online way and then improve the well-known PETRELS algorithm to update the underlying subspace in the missing data context.


# DEMO

Run the file DEMO_SEP_Main.m for similated data

Run the file DEMO_Video.m for real video dataset

The Lobby data: https://drive.google.com/drive/folders/11a_TgkJAyw7PvF-lz9RuUW_SHeMk_F1H?usp=sharing 

# State-of-the-art algorithms for comparison
+ GRASTA: https://sites.google.com/site/hejunzz/grasta
+ ROSETA: http://www.merl.com/research/license#ROSETA
+ ReProCS: https://github.com/praneethmurthy/ReProCS
+ NORST: https://github.com/praneethmurthy/NORST

# References

This code is free and open source for research purposes.
If you use this code, please acknowledge the following papers.

[1] L.T. Thanh, V.D. Nguyen, N. L. Trung and K. Abed-Meraim. “Robust Subspace Tracking with Missing Data and Outliers: Novel Algorithm with Convergence Guarantee”. IEEE
Trans. Signal Process., 2021. (to appear)

[2] L.T. Thanh, V.D Nguyen, N.L. Trung and K. Abed-Meraim. “Robust Subspace Tracking with Missing Data and Outliers via ADMM”. European Signal Process. Conf. (EUSIPCO), 2019


