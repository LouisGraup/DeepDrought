# use Segment Anything model to create masks for hemispherical photos

import numpy as np
import matplotlib.pyplot as plt
import torch
import cv2

# helper functions
def show_mask(mask, ax, random_color=False):
    if random_color:
        color = np.concatenate([np.random.random(3), np.array([0.6])], axis=0)
    else:
        color = np.array([30/255, 144/255, 255/255, 0.6])
    h, w = mask.shape[-2:]
    mask_image = mask.reshape(h, w, 1) * color.reshape(1, 1, -1)
    ax.imshow(mask_image)

def show_points(coords, labels, ax, marker_size=375):
    pos_points = coords[labels==1]
    neg_points = coords[labels==0]
    ax.scatter(pos_points[:, 0], pos_points[:, 1], color='green', marker='*', s=marker_size, edgecolor='white', linewidth=1.25)
    ax.scatter(neg_points[:, 0], neg_points[:, 1], color='red', marker='*', s=marker_size, edgecolor='white', linewidth=1.25)   
    
def show_box(box, ax):
    x0, y0 = box[0], box[1]
    w, h = box[2] - box[0], box[3] - box[1]
    ax.add_patch(plt.Rectangle((x0, y0), w, h, edgecolor='green', facecolor=(0,0,0,0), lw=2))    


image = cv2.imread('ISF2020/2021-11.png')

#plt.figure(figsize=(10,10))
plt.imshow(image)
plt.axis('on')
plt.show()


from segment_anything import sam_model_registry, SamPredictor
sam = sam_model_registry["vit_h"](checkpoint="sam_vit_h_4b8939.pth")
#sam.to(device="cuda")

predictor = SamPredictor(sam)
predictor.set_image(image)


input_point = np.array([[1400, 1050]])
input_label = np.array([1])

plt.imshow(image)
show_points(input_point, input_label, plt.gca())
plt.axis('on')
plt.show()  


masks, scores, logits = predictor.predict(
    point_coords=input_point,
    point_labels=input_label,
    multimask_output=True,
)


for i, (mask, score) in enumerate(zip(masks, scores)):
    plt.figure(figsize=(10,10))
    plt.imshow(image)
    show_mask(mask, plt.gca())
    show_points(input_point, input_label, plt.gca())
    plt.title(f"Mask {i+1}, Score: {score:.3f}", fontsize=18)
    plt.axis('off')
    plt.show()  

# Save the first mask as an example
mask = masks[1]
mask = mask.astype(np.uint8) * 255
#mask[mask==0] = np.nan
cv2.imwrite('mask20_63.png', mask)