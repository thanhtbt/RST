# PETRELS-ADMM: Robust Subspace Tracking with Missing Data and Outliers

We propose a novel algorithm called PETRELS-ADMM to deal with subspace tracking in the presence of outliers and missing data. The proposed approach consists of two main stages: outlier rejection and subspace estimation. Particularly, we first use ADMM solver for detecting outliers living in the measurement data in an efficient online way and then improve the well-known PETRELS algorithm to update the underlying subspace in the missing data context.

## Updates:
+ Jan 2021: Create this repository.

+ Oct 16th 2021: Reorganize the entire repository.

## DEMO

+ Run "DEMO_SEP_Main.m" for synthetic data.

+ Run "DEMO_Video.m" for real data: The Lobby video data can be downloaded from the Release or [here](https://drive.google.com/drive/folders/11a_TgkJAyw7PvF-lz9RuUW_SHeMk_F1H?usp=sharing) or Release.

## State-of-the-art algorithms for comparison
+ GRASTA: https://sites.google.com/site/hejunzz/grasta
+ ROSETA: http://www.merl.com/research/license#ROSETA
+ ReProCS: https://github.com/praneethmurthy/ReProCS
+ NORST: https://github.com/praneethmurthy/NORST

## Some results

Similated data: matrix completion and performance comparsion between PETRELS-ADMM and the state-of-the-art RST algorithms
<p float="left">
  <img src="https://user-images.githubusercontent.com/26319211/110497389-b01e0680-80f6-11eb-870f-88e9b6c65ae4.PNG" width="300" height='250' />
  <img src="https://user-images.githubusercontent.com/26319211/110496361-af38a500-80f5-11eb-84c3-26485bda807c.jpg" width="300" height='250' /> 
</p>

Video background-foreground separation application

<img src="https://user-images.githubusercontent.com/26319211/110496363-afd13b80-80f5-11eb-8770-30510d66e271.PNG" width="300" height='250'>


## References

This code is free and open source for research purposes. If you use this code, please acknowledge the following papers.

[1] L.T. Thanh, V.D. Nguyen, N. L. Trung and K. Abed-Meraim. “[*Robust Subspace Tracking with Missing Data and Outliers: Novel Algorithm with Convergence Guarantee*](https://drive.google.com/file/d/1MIFmZlQyx1L3lUZ5QGz6ARWVPD4eBYgF/view?usp=sharing)”. **IEEE Trans. Signal Process., 69:2070–2085, 2021**. [[DOI](https://ieeexplore.ieee.org/document/9381678)],[[PDF](https://drive.google.com/file/d/1MIFmZlQyx1L3lUZ5QGz6ARWVPD4eBYgF/view?usp=sharing)].

[2] L.T. Thanh, V.D Nguyen, N.L. Trung and K. Abed-Meraim. “[*Robust Subspace Tracking with Missing Data and Outliers via ADMM*](https://drive.google.com/file/d/1fOfWjUdMgUuOI7yWpouid3BMb29QQzkr/view?usp=sharing)”. **27th European Signal Process. Conf. (EUSIPCO), 1-5,2019**. [[DOI]( https://ieeexplore.ieee.org/document/8903031)],[[PDF](https://drive.google.com/file/d/1fOfWjUdMgUuOI7yWpouid3BMb29QQzkr/view?usp=sharing)].

