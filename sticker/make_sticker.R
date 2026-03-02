# File to make and adjust the hexagonal sticker

# load the necessary packages
library(hexSticker) # hexSticker generator
library(magick)     # Advanced image processing
library(sysfonts)   # font selection
library(tidyverse)

#setwd("/Users/tijnjacobs/Library/CloudStorage/OneDrive-VrijeUniversiteitAmsterdam/Documents/GitHub/ShrinkageTrees/sticker")

pine_img <- magick::image_read('vector_graph_logo_high.png')

fonts_df <- font_files()

font_add("Avenir Next", "Avenir Next.ttc")
font_add("Avenir Next Condensed", "Avenir Next Condensed.ttc")
font_add("DIN Condensed", "DIN Condensed Bold.ttf")
font_add("DIN Alternate", "DIN Alternate Bold.ttf")
font_add("Helvetica Neue", "HelveticaNeue.ttc")


sticker(
  
  # title of the package
  package = "ShrinkageTrees",
  p_size = 30,
  p_y = 0.55,
  
  # Border thickness
  h_size = 1.4,
  
  # Image
  subplot = pine_img,  
  s_x = 1,
  s_y = 1.2,
  s_width = 1.35,
  s_height = 1.35,
  
  # Colouring
  p_color = "#3F7426",
  h_fill = "#B9DA69",  # background color (feel free to change)
  h_color = "#3F7426", # border color
  
  # Spotlight
  spotlight = TRUE,
  l_x = 1.35,
  l_y = 1.55,
  l_width = 2,
  l_alpha = 0.21,
  
  filename = "ShrinkageTrees_hex.png",
  dpi = 600,
  p_family = "DIN Alternate"
)
