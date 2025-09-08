These codes are intended to train a UNET model + three different configuration of GAN-based models for downscaling CMIP6 projections from ~200 km to ~50 km. 

Model 1. The UNET model (UNET)  served as a baseline for comparison. It employs an encoder-decoder structure with skip connections to capture spatial features of precipitation fields. The encoder reduces spatial dimensions while extracting hierarchical features, and the decoder reconstructs high-resolution outputs. This architecture is widely used for image-to-image translation tasks, and recently, it dominated the field of weather forecasting (e.g., Kaparakis and Mehrkanoon, 2023; Li et al., 2024; Trebing et al., 2021).

Model 2. The second model (WGAN)  integrates a ConvLSTM layer into the generator of the WGAN. The generator maps low-resolution precipitation inputs to high-resolution outputs while capturing temporal dependencies. The discriminator critic evaluates the quality of generated fields without temporal layers that explicitly understand the temporal dependence. However, the discriminator critic uses a set of fully connected layers that decompose and process the third dimension (i.e., time).

Model 3. To address the challenge of intermittent precipitation, the third model (WGAN_NZ) includes a nonzero thresholding layer within the generator in addition to the structure of the previous model. This layer ensures that generated precipitation fields do not include negative values, even or positive ones very closer to zero, effectively distinguishing between zero and nonzero precipitation regions. 

Model 4. The fourth model  (WGAN_NZ_DT) combines the advantages of temporal and threshold-aware mechanisms. ConvLSTM layers were integrated into both the generator and discriminator the critic to capture spatiotemporal dependencies. The nonzero thresholding layer in the generator further enforces the realism of sparse precipitation patterns.



The trained models are eventually saved in the "Outputs" directory.

The evaluation codes (in R) are saved in "EvaluationCodes" directory.




Data Availability: The inputs for the training codes are the CanRCM5 (Canadian Regional Climate model) historical simulation. The format of these file should be (Latitude, Longitude, hours) where each input file has 24 time steps representing one day entery of a field. The Canadian Praries is the study region (The geographic extent of the study ranges from 121° West to 94° West and from 48° North to 61° North). To downscale CMIP6 projects, the simulations should be reshaped to the same format.

Data Usage:
There are several ways to utilize the provided codes and models:

1. Use Pre-Trained Models: The trained models in the "Outputs" folder can be used to downscale any climate model output from 2° × 2° (~200 km) to 0.5° × 0.5° (~50 km). Note that these models were trained using the Canadian Regional Model (CanRCM4), and while transfer learning is not guaranteed, it may still be worth exploring.

2. Fine-Tune for New Data: The trained models in the "Outputs" folder can serve as pre-trained models, potentially enabling faster convergence when trained on new datasets.

3. Train from Scratch: Users can select any preferred model from the code, modify the data directory to their target dataset, and train the model from scratch.

Additional Notes:
* The models are trained using paired low- and high-resolution data. Once trained, they take low-resolution fields as input and generate downscaled, high-resolution fields.
* The structure of the low- and high-resolution fields is predefined within the models:
    - Low-resolution fields: 14 longitudes × 7 latitudes
    - High-resolution fields: 54 longitudes × 26 latitudes
    - If the structure differs for a user’s dataset, the code must be adjusted accordingly.
* The models incorporate a temporal component, meaning they are designed to process 24 sequential fields at a time.