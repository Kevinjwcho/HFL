# Estimating Historical Functional Linear Models with a Nested Group Bridge Approach

This repository consists of coding for simulation studies and application in the manuscript entitled "Estimating Historical Functional Linear Models with a Nested Group Bridge Approach".

The folder "simulation" consists of:

1. Folder "datafunction": include all the functions used to generate simulation data, realize our proposed nested group bridge method and the truncation methods A and B;

2. Folder "beta1": include R code for the simulation study for scenario I in the manuscript. The folder includes simulation results for the proposed method, the smoothing 
         spline mehtod and the truncation methods A and B. 

         The R codes for the proposed method and the smoothing spline method are in the folder "NGR": 
           (i)   bt1_NGR_n100_SNR2.R: the R code for the proposed method with the simulation setting of scenario I, signal-to-noise ratio 2, and sample size 100;
           (ii)  bt1_NGR_n500_SNR2.R: the R code for the proposed method with the simulation setting of scenario I, signal-to-noise ratio 2, and sample size 500;
           (iii) bt1_ss_n100_SNR2.R: the R code for the smoothing spline method with the simulation setting of scenario I, signal-to-noise ratio 2, and sample size 100;
           (iv)  bt1_ss_n100_SNR2.R: the R code for the smoothing spline method with the simulation setting of scenario I, signal-to-noise ratio 2, and sample size 500;

         The R codes for the truncation method are in the folder "tr": 
           (i)  bt1_n100_SNR2: the R code for the truncation methods with the simulation setting of scenario I, signal-to-noise ratio 2, and sample size 100;
           (ii) bt1_n500_SNR2: the R code for the truncation methods with the simulation setting of scenario I, signal-to-noise ratio 2, and sample size 500;
          
3. Folder "beta2": include R code for the simulation study for scenario II in the manuscript. The R codes are exactly the same as the ones in folder "beta1", except that the codes in 
         this folder are for scenario II.

4. Folder "beta3": include R code for the simulation study for scenario III in the manuscript. The R codes are exactly the same as the ones in folder "beta1", except that the codes in 
         this folder are for scenario III.



The folder "application" consists of:
1.Truck1Run1.csv: the truck particulate emissions data
2.truck_analysis.R: the R code for analyzing the particulate emissions data. We construct the bootstrap confidence intervals by resampling the residuals. We also compare oru proposed 
                    method with the smoothing spline method.
