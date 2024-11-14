#!/bin/bash

# Create necessary folders
mkdir -p cropped_images
mkdir -p superposed

echo "Starting image processing..."

# Step 1: Copy and crop each image, saving the cropped version in "cropped_images"
echo "Copying and cropping images..."
for img in Figure_*.png; do
    cp "$img" "cropped_images/$img"  # Copy to cropped_images folder
    mogrify -gravity Center -crop 50%x40%+120-100 "cropped_images/$img"  # Crop the copy
done

echo "Image cropping completed. Proceeding to superpose pairs..."

# Step 2: Superpose each pair of cropped images
i=1
while true; do
    img1="cropped_images/Figure_${i}.png"
    img2="cropped_images/Figure_$((i+1)).png"
    
    # Check if both images in the pair exist
    if [[ -f "$img1" && -f "$img2" ]]; then
        output="superposed/Figure_${i}_superposed.png"
        
        echo "Superposing $img1 and $img2 into $output"
        
        # Superpose the images with 50% transparency
        composite -dissolve 50% "$img1" "$img2" "$output"
        
        # Check if the output image was created successfully
        if [[ -f "$output" ]]; then
            echo "Successfully created $output"
        else
            echo "Failed to create $output"
        fi
        
        # Increment by 2 to move to the next pair
        i=$((i+2))
    else
        # Exit loop if we run out of pairs
        echo "No more pairs to process. Exiting loop."
        break
    fi
done

echo "Image superposition completed. Superposed images saved in the 'superposed' folder."
