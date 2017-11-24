convert -crop 4x2000+0+0 mode.pbm mode_crop.pbm
mogrify -filter Point -resize 40x2000 mode_crop.pbm
convert -crop 1403x2000+0+0  mask.pbm mask_crop.pbm
convert -crop 1403x2000+0+0  dist.pbm dist_crop.pbm
montage -tile 3x1 -geometry +10+10 mode_crop.pbm dist_crop.pbm mask_crop.pbm  montage.png
