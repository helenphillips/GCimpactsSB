# dat$TaxaGroupGSBA <- 
  
dat$TaxaGroup[grepl("Arthropods \\(mites,", ignore.case = TRUE, dat$TaxaGroup)] <- "Soil arthropods"

# microfauna

grepl("Nematode", ignore.case = TRUE, dat$TaxaGroup), "Nematodes")

grepl("Rotifer", ignore.case = TRUE, dat$TaxaGroup), "Rotifers")

grepl("Tardigrade", ignore.case = TRUE, dat$TaxaGroup), "Tardigrades")


# mesofauna

grepl("Protura", ignore.case = TRUE, dat$TaxaGroup), "Protura")

grepl("Collembola|springtail", ignore.case = TRUE, dat$TaxaGroup), "Collembola")

grepl("Enchytraeid", ignore.case = TRUE, dat$TaxaGroup), "Enchytraeids")

grepl("Acari|stigmat|Oribati|Total mites|predatory mites|Mites", ignore.case = TRUE, dat$TaxaGroup), "Acari")

grepl("Diplur", ignore.case = TRUE, dat$TaxaGroup), "Diplura")

grepl("Pseudoscorpio", ignore.case = TRUE, dat$TaxaGroup), "Pseudoscorpions")


# macrofauna

grepl("ants|formicid", ignore.case = TRUE, dat$TaxaGroup), "Formicidae")

grepl("termites", ignore.case = TRUE, dat$TaxaGroup), "Termites")

grepl("isopod", ignore.case = TRUE, dat$TaxaGroup), "Isopoda")

grepl("milliped|myriapod|centiped|pauropod|symphyl", ignore.case = TRUE, dat$TaxaGroup), "Myriapoda")

grepl("earthworm", ignore.case = TRUE, dat$TaxaGroup), "Earthworms")

grepl("beetle|coleopter", ignore.case = TRUE, dat$TaxaGroup), "Coleoptera")

grepl("Larvae", ignore.case = TRUE, dat$TaxaGroup), "Soil Insect Larvae")


# ground or litter-dwellers macrofauna

grepl("spider|arachnid|aranea", ignore.case = TRUE, dat$TaxaGroup), "Spiders")

grepl("snail|slug|gastropod", ignore.case = TRUE, dat$TaxaGroup), "Gastropoda")

grepl("bee|wasp", ignore.case = TRUE, dat$TaxaGroup), "Ground-nesting bees")


grepl("microarthropod|micro-arthropod|meso-arthropod", ignore.case = TRUE, dat$TaxaGroup), "Micro-arthropods")


grepl("macro-arthrop|macroarthrop", ignore.case = TRUE, dat$TaxaGroup), "Macro-arthropods")



