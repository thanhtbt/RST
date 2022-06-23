# PETRELS-ADMM: Robust Subspace Tracking with Missing Data and Outliers

We propose a novel algorithm called PETRELS-ADMM to deal with subspace tracking in the presence of outliers and missing data. The proposed approach consists of two main stages: outlier rejection and subspace estimation. Particularly, we first use ADMM solver for detecting outliers living in the measurement data in an efficient online way and then improve the well-known PETRELS algorithm to update the underlying subspace in the missing data context.

## Updates:
+ Jan 2021: Create this repository.

+ Oct 2021: Reorganize the entire repository.

+ Jun 2022: Add a demo DEMO_SEP_Main_Synthetic.m (avoids the warning message from MATLAB 202x caused by ReProCS and NORST)

## DEMO

+ Run "DEMO_SEP_Main.m" (Mathlab R201x) or "DEMO_SEP_Main_Synthetic.m" (Mathlab R202x) for synthetic data .

+ Run "DEMO_Video.m" for real data: Video data can be downloaded from the Release or [here](https://drive.google.com/drive/folders/11a_TgkJAyw7PvF-lz9RuUW_SHeMk_F1H?usp=sharing).

## State-of-the-art algorithms for comparison
+ GRASTA: https://sites.google.com/site/hejunzz/grasta
+ ROSETA: http://www.merl.com/research/license#ROSETA
+ ReProCS: https://github.com/praneethmurthy/ReProCS
+ NORST: https://github.com/praneethmurthy/NORST

## Some results

Similated data:  

![untitled](https://user-images.githubusercontent.com/26319211/175233846-5d1a564a-9057-498c-9b10-66d4469578e6.jpg)

Video background-foreground separation application

<img src="https://user-images.githubusercontent.com/26319211/110496363-afd13b80-80f5-11eb-8770-30510d66e271.PNG" width="700" height='550'>


## References

This code is free and open source for research purposes. If you use this code, please acknowledge the following papers.

[1] L.T. Thanh, V.D. Nguyen, N. L. Trung and K. Abed-Meraim. “[*Robust Subspace Tracking with Missing Data and Outliers: Novel Algorithm with Convergence Guarantee*](https://drive.google.com/file/d/1MIFmZlQyx1L3lUZ5QGz6ARWVPD4eBYgF/view?usp=sharing)”. **IEEE Trans. Signal Process., 69:2070–2085, 2021**. [[DOI](https://ieeexplore.ieee.org/document/9381678)],[[PDF](https://drive.google.com/file/d/1MIFmZlQyx1L3lUZ5QGz6ARWVPD4eBYgF/view?usp=sharing)].

[2] L.T. Thanh, V.D Nguyen, N.L. Trung and K. Abed-Meraim. “[*Robust Subspace Tracking with Missing Data and Outliers via ADMM*](https://drive.google.com/file/d/1fOfWjUdMgUuOI7yWpouid3BMb29QQzkr/view?usp=sharing)”. **27th European Signal Process. Conf. (EUSIPCO), 1-5,2019**. [[DOI]( https://ieeexplore.ieee.org/document/8903031)],[[PDF](https://drive.google.com/file/d/1fOfWjUdMgUuOI7yWpouid3BMb29QQzkr/view?usp=sharing)].

