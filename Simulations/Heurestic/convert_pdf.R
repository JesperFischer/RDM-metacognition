
library(magick)
# Path to your PDF
pdf_path <- here::here("Simulations","Heurestic","pdfs")



pdf_files <- list.files(pdf_path, pattern = "\\.pdf$", full.names = TRUE)
pdf_names <- list.files(pdf_path, pattern = "\\.pdf$")

i = 0

for(pdf in pdf_files){
  
  i = i +1
  # Read PDF at 300 DPI
  img <- magick::image_read_pdf(pdf, density = 300)
  
  # Write to TIFF (one file per page if multi-page PDF)
  magick::image_write(img, path = here::here("Simulations","Heurestic","tiffs",paste0(pdf_names[i],".tiff")), format = "tiff")
  
}
