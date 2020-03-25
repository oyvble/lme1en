# lme1en
A linear model with random intercept effects and Elastic Net penalty 

Suitable for (p >> n)


Installation:
devtools::install_github("oyvble/lme1en")


Modell: 
For observation i within batch b:

![y_(b,i)=x_(b,i)^T β+η_b+ϵ_(b,i)](https://render.githubusercontent.com/render/math?math=y_%7Bb%2C%20i%7D%3D%5Cboldsymbol%7Bx%7D_%7Bb%2C%20i%7D%5E%7BT%7D%20%5Cboldsymbol%7B%5Cbeta%7D%2B%5Ceta_%7Bb%7D%2B%5Cepsilon_%7Bb%2C%20i%7D)

![\epsilon_{b, i} \sim \text {iid. } N\left(0, \sigma_{\epsilon}^{2}\right)](https://render.githubusercontent.com/render/math?math=%5Cepsilon_%7Bb%2C%20i%7D%20%5Csim%20%5Ctext%20%7Biid.%20%7D%20N%5Cleft%280%2C%20%5Csigma_%7B%5Cepsilon%7D%5E%7B2%7D%5Cright%29%0A)

![\eta_{b} \sim \text {iid.} N\left(0, \sigma_{\eta}^{2}\right)](https://render.githubusercontent.com/render/math?math=%5Ceta_%7Bb%7D%20%5Csim%20%5Ctext%20%7Biid.%7D%20N%5Cleft%280%2C%20%5Csigma_%7B%5Ceta%7D%5E%7B2%7D%5Cright%29%0A)

Elastic Net penalties on beta. The package use following reparameterisation:

![\rho = \frac{\sigma_{\eta}^{2}}{\sigma_{\epsilon}^{2}+\sigma_{\eta}^{2}}](https://render.githubusercontent.com/render/math?math=%5Crho%20%3D%20%5Cfrac%7B%5Csigma_%7B%5Ceta%7D%5E%7B2%7D%7D%7B%5Csigma_%7B%5Cepsilon%7D%5E%7B2%7D%2B%5Csigma_%7B%5Ceta%7D%5E%7B2%7D%7D%0A)
