#!/bin/bash

montage -geometry +1+1 x/* x_out.png
montage -geometry +1+1 y/* y_out.png
montage -geometry +1+1 z/* z_out.png
echo "Images joined!"