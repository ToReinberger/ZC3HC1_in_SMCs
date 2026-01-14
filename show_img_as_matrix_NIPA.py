import os

import PIL.Image
import matplotlib.pyplot as plt
import glob
import PIL
import numpy as np
import pandas as pd
from tkinter import filedialog


# layout = fr"Z:\Projects\Neointima_ZC3HC1_Cells\Neointima_formation\Zuordnung_NIPA_Femuralis.xlsx"
# allocation = pd.read_excel(layout)
# print(allocation)
path = r"20230228_mSMC_ZC3_KO_Ccbn1\*\*cropped*.png"
# tag = "mSMC_Ccbn1_NIPA_KO_NEW"
# tag = "mSMC_Ccbn1_NIPA_KO_NEW"
tag = "mSMC_CCNB1_NIPA_KO_NEW"
files = glob.glob(path)
images = []

dim2_cols = int(len(files) / 4)
"""dim1_rows = int(len(images) / dim2_cols)
if dim2_cols * dim1_rows < len(images):
    dim1_rows += 1"""
dim1_rows = 4
if dim1_rows * dim2_cols < len(images):
    print("increase space")
    dim1_rows = dim1_rows * 2

print(dim1_rows, dim2_cols)
dim1_row = 0
dim2_col = 0
data = []
genotype = ""
for img_path in files:
    print(img_path)
    if "TEMP" in img_path:
        continue
    # print(mouse_id, genotype, treatment)
    # images.append(img_path.split("\\")[-1])
    image = PIL.Image.open(img_path)
    images.append(image)
    """x, y = image.size
    image = np.asarray(image)
    crop = 0.15  # 10%
    image = image[int(crop * y):int((1 - crop) * y), int(crop * x):int((1 - crop) * x)]"""

# tbl = np.asarray(tbl)
# tbl = sorted(tbl, key=lambda x: x[0])

# plt.figure(figsize=(12, 18),# layout="tight"
#           )
# plt.subplots_adjust(left=0.025, right=0.975, bottom=0.1, top=0.95, wspace=0.01, hspace=0.01)
plt.figure(figsize=(dim2_cols * 3, 12))
plt.subplots_adjust(left=0.02, right=0.98, top=.98, bottom=0.02, wspace=0.03, hspace=0.03)

for image in images:

    plt.subplot2grid((dim1_rows, dim2_cols), (dim1_row, dim2_col))
    # image = np.asarray(image)
    plt.imshow(image)
    plt.axis("off")
    # plt.title(genotype + "," + treatment, fontdict={"size": 8})
    # plt.show()
    dim1_row += 1
    # dim1_row += 1
    if dim1_row == dim1_rows:
        dim1_row = 0
        dim2_col += 1


plt.savefig(rf"C:\Users\tobia\PycharmProjects\WetLab\Image_Analyzer\Confocal_IF\{tag}.pdf")
plt.savefig(rf"C:\Users\tobia\PycharmProjects\WetLab\Image_Analyzer\Confocal_IF\{tag}.tif")
os.startfile(rf"C:\Users\tobia\PycharmProjects\WetLab\Image_Analyzer\Confocal_IF\{tag}.tif")
# plt.show()
quit()
