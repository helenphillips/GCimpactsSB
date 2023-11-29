# This is a ridiculously long script just to make one supplementary figure
# Uses the model of taxa group X stressor.
# As there was one factor level missing (earthworms in flooding), some pre-made functions
# I found did not work
# in hindsight, probably should have just done four separate models.



library(Hmisc)
library(viridis)

mod.gsba.taxa_stressor <- readRDS(file = "Models/GSBAMod_stressors.rds")


# summary(mod.gsba.taxa_stressor)

newdat <- predict(mod.gsba.taxa_stressor, 
        
        rbind( # intercept (compost)
          
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            1,0,0, # collembola
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,1,0, # earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),      
          
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),   
          
          
          # fire
          c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),      

          c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            1,0,0, # collembola
            1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          
          c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,1,0, # earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),      
          
          c(1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),      
          
          # co2
          c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),      
          
          c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            1,0,0, # collembola
            0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          
          c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,1,0, # earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),      
          
          c(0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0),   
          
          
          # now just for acari
          # grazing
          c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
          # harvesting
          c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
          # Manure + Slurry
          c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
          # Metals
          c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Organic versus Inorganic
          c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Other Organic
          c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Pesticides
          c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Residue + Mulch
          c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Sludge
          c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Synthetic Fertilizers 
          c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # STemperature
          c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Tillage
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # drought
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # flood
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
            0,0,0, # acari
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          
          
          # now just for collembola
          # grazing
          c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
            1,0,0, # collembola
            0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
          # harvesting
          c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
            1,0,0, # collembola
            0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
          # Manure + Slurry
          c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
            1,0,0, # collembola
            0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
          # Metals
          c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
            1,0,0, # collembola
            0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Organic versus Inorganic
          c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
            1,0,0, # collembola
            0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Other Organic
          c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
            1,0,0, # collembola
            0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Pesticides
          c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
            1,0,0, # collembola
            0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Residue + Mulch
          c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
            1,0,0, # collembola
            0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Sludge
          c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
            1,0,0, # collembola
            0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Synthetic Fertilizers 
          c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
            1,0,0, # collembola
            0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # STemperature
          c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
            1,0,0, # collembola
            0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Tillage
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
            1,0,0, # collembola
            0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # drought
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
            1,0,0, # collembola
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # flood
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
            1,0,0, # collembola
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          
          
          
          # now just for Earthworms
          # grazing
          c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,1,0, # Earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
          # harvesting
          c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
            0,1,0, # Earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
          # Manure + Slurry
          c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
            0,1,0, # Earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0), 
          # Metals
          c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
            0,1,0, # Earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Organic versus Inorganic
          c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
            0,1,0, # Earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Other Organic
          c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
            0,1,0, # Earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Pesticides
          c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
            0,1,0, # Earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Residue + Mulch
          c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
            0,1,0, # Earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Sludge
          c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
            0,1,0, # Earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Synthetic Fertilizers 
          c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
            0,1,0, # Earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # STemperature
          c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
            0,1,0, # Earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # Tillage
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
            0,1,0, # Earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # drought
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
            0,1,0, # Earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          # flood
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
            0,1,0, # Earthworms
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
          
          
          # now just for nematodes
          # grazing
          c(0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0), 
          # harvesting
          c(0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0), 
          # Manure + Slurry
          c(0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0), 
          # Metals
          c(0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0),
          # Organic versus Inorganic
          c(0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0),
          # Other Organic
          c(0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0),
          # Pesticides
          c(0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0),
          # Residue + Mulch
          c(0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0),
          # Sludge
          c(0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0),
          # Synthetic Fertilizers 
          c(0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0),
          # STemperature
          c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0),
          # Tillage
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0),
          # drought
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0),
          # flood
          c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,
            0,0,1, # nematodes
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1)
        ),
        addx=TRUE, digits=2
)


gcdlevels <- c(rep("compost", times = 4), 
               rep("fire", times = 4),  
               rep("gas co2", times = 4),
rep(c("Grazing","Harvesting" , "Manure + Slurry"    ,                                       
"Metals","Organic versus Inorganic",
"Other Organic fertilisers (NOT including compost and Urea)",
"Pesticides","Residue + Mulch","Sludge (including Biosolids)",                         
"Synthetic Fertilizers" ,"Temperature"  ,"Tillage",                                             
"WaterAvailability-Drought", "WaterAvailability-Flood"), times = 4))



taxa <- c(rep(c("acari", "collembola", "earthworms", "nematodes"), times = 3),
          rep(c("acari", "collembola", "earthworms", "nematodes"), each = 14))


## order
order <- c(
# flooding
which(gcdlevels == "WaterAvailability-Flood"),
# drought
which(gcdlevels == "WaterAvailability-Drought"),
# temp
which(gcdlevels == "Temperature"),
#co2
which(gcdlevels == "gas co2"),

# grazing
which(gcdlevels == "Grazing"),
# fire
which(gcdlevels == "fire"),
# harvesting
which(gcdlevels == "Harvesting"),
# conventional farming
which(gcdlevels == "Organic versus Inorganic"),
# tillage
which(gcdlevels == "Tillage"),


# pesticides
which(gcdlevels == "Pesticides"),
# metals
which(gcdlevels == "Metals"),

# synthetic
which(gcdlevels == "Synthetic Fertilizers"),
# compost
which(gcdlevels == "compost"),
# sludge
which(gcdlevels == "Sludge (including Biosolids)"),
# manure
which(gcdlevels == "Manure + Slurry"),
# residue
which(gcdlevels == "Residue + Mulch"),
#other
which(gcdlevels == "Other Organic fertilisers (NOT including compost and Urea)")
)


oorder <- rev(order)
gcdlevels[which(gcdlevels == "Other Organic fertilisers (NOT including compost and Urea)")] <- "Other Organic \nfertilisers"

gcdlevels[which(gcdlevels == "WaterAvailability-Flood")] <- "Flooding"
gcdlevels[which(gcdlevels == "WaterAvailability-Drought")] <- "Drought"
gcdlevels[which(gcdlevels == "gas co2")] <- "Carbon dioxide \nincrease"
gcdlevels[which(gcdlevels == "Temperature")] <- "Temperature \nchange"
gcdlevels[which(gcdlevels == "fire")] <- "Fire"
gcdlevels[which(gcdlevels == "Organic versus Inorganic")] <- "Conventional \nfarming"
gcdlevels[which(gcdlevels == "compost")] <- "Compost"
gcdlevels[which(gcdlevels == "Sludge (including Biosolids)")] <- "Sludge"

gcdlevels[which(gcdlevels == "Residue + Mulch")] <- "Residue + \nMulch"
gcdlevels[which(gcdlevels == "Synthetic Fertilizers")] <- "Synthetic \nfertilizers"
gcdlevels[which(gcdlevels == "Manure + Slurry")] <- "Manure + \nSlurry"



stressors <- gcdlevels[oorder]
stressors <- stressors[seq(1, 68, by = 4)]

# put the data into order now
newdat <- newdat[oorder]

# fake df to get a row index (because I just can't figure out the str of the other one)
testdf <- as.data.frame(newdat$X)
# this row technically shouldn't exist, as no underlying data
indx <- which(testdf[,'GCDTypeWaterAvailability-Flood'] == 1 & testdf[,'GSBAEarthworms'] == 1)

newdat <- newdat[-indx,]



taxalabs <- rev(rep(c("Acari", "Collembola", "Earthworms", "Nematodes"), times = 17))
taxalabs <- taxalabs[-indx]
# This is because I can't get the left margin to be bigger...annoying
taxalabs <- paste("                          ", taxalabs)
# labs <- paste(gcdlevels[oorder],taxalabs)





cols <- rep(viridis(4), length = 68) # each taxa has a colour
cols <- cols[-indx]

cols[which(newdat$ci.ub > 0 & newdat$ci.lb < 0)] <- "#808080FF" # not significant

# "#d8d8d8FF"

pdf(file = file.path("C:/Users/helenp/WORK/GCimpactsSB/Figures", "Taxaxstressors.pdf"),
    width = 12, height = 24)

par(mar=c(5, 10, 2, 8))


errbar(x =taxalabs, y = newdat$pred, 
       yplus = newdat$ci.ub,
       yminus=newdat$ci.lb, 
       cex = 2,
       col = cols,
       xaxt="n")


polygon(x=c(-3,-3,3,3), y=c(24.5,32.5,32.5,24.5), col="#DCDCDC7D", border=F)
polygon(x=c(-3,-3,3,3), y=c(52.5,67.5,67.5,52.5), col="#DCDCDC7D", border=F)


# repeat over the shading
par(new = TRUE)

errbar(x =taxalabs, y = newdat$pred, 
       yplus = newdat$ci.ub,
       yminus=newdat$ci.lb,
       cex = 2,
       col = cols,
       xaxt="n")



abline(v=0, lty = 2)

ats <- (seq(2.5, 67, by = 4))
# manually override last one, as we are missing a category
ats[17] <- 66
for(a in 1:length(stressors)){
  mtext(stressors[a], side = 2, line = 4.5, at = ats[a], cex = 1)
}




mtext("Climate change", side = 2, line = 7.5, at = 59.5, cex = 1.2)
mtext("Land-use intensification", side = 2, line = 7.5, at = 42.5, cex = 1.2)
mtext("Pollution", side = 2, line = 7.5, at = 28.5, cex = 1.2)
mtext("Nutrient enrichment", side = 2, line = 7.5, at = 12.5, cex = 1.2)



mtext("Effect size", side = 1, line = 3, cex = 1.5)


## The n's papers
tapply(taxa_dat2$ID, list( taxa_dat2$GSBA, taxa_dat2$GCDType), function(x) length(unique(x)))
table(taxa_dat2$GSBA, taxa_dat2$GCDType)


caseNs <- c(12, 16, 6, 10, # other organic
            18, 21, 15, 21, # residue
            27, 36, 7, 23, # manure
            10, 3, 9, 39, # sludge 
            12, 17, 7, 18, # compost
            58, 35, 38, 55, # synthetic
            122, 35, 39, 37, # metals
            70, 44, 44, 50, # pesticides
            44, 76, 21, 31 , # tillage
            39, 27, 4, 13, # inorganic
            24, 15, 16, 21, # harvesting
            12, 4, 19, 37, # fire
            49, 28, 8, 9, # grazing
            29, 11, 15, 23, # co2
            27, 4, 36, 72, # temperature
            16, 4, 15, 39, # drought
            23, 8, 7 )# flood (missing earthworms)
            
paperNs <- c(5, 4, 4, 4, # other organic
             12, 9, 19, 9, # residue
             20, 16, 3, 5, # manure
             5, 2, 4, 6, # sludge
             7, 3, 1, 2, # compost
             37, 19, 23, 21, # synthetic
             38, 15, 19, 14, # metals
             18, 12, 14, 14, # pesticides
             29, 34, 12, 12, # tillage
             20, 15, 4, 5, # inorganic
             12, 6, 9, 9, # harvesting
             4, 3, 12, 14, # fire
             20, 14, 6, 5, # grazing
             14, 2, 10, 8, # co2
             15, 3, 20, 17, # temperature
             10, 2, 10, 13, #drought
             8, 6, 4) # flood (missing earthworms)

for(n in 1:length(paperNs)){
  mtext(paste0("n = ", caseNs[n], " (", paperNs[n], ")"), side = 4, line = 0, at = n + 0.1, las = 2)
}

# remove abline from top
polygon(x=c(-0.5,-0.5,0.4,0.4), y=c(67.5,72,72,67.5), col="white", border=F)



legendcols <- c(rev(viridis(4)),"#808080FF")

legend("topright", inset=c(0.002,0.002), legend=c("Acari", "Collembola", "Earthworms", "Nematodes", "Not sign."), pch = 19,
       horiz=TRUE, col = legendcols, cex = 1, bty = "o", bg = "white", pt.cex = 1.5)






dev.off()