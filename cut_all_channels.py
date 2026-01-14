import PIL
import glob
import matplotlib.pyplot as plt
from PIL import Image, ImageEnhance


path = r"C:\Users\tobia\PycharmProjects\WetLab\Image_Analyzer\Confocal_IF"
file = "20230228_mSMC_ZC3_KO_Ccbn1"
phase = "10_Migr"

for img_path in glob.glob(path + rf"\{file}\{phase}\*super_projection_CH*.png"):
    if "cropped" in img_path:
        continue
    print(img_path)
    img = Image.open(img_path)

    img = img.rotate(11, PIL.Image.Resampling.NEAREST, expand=0)
    img = ImageEnhance.Brightness(img).enhance(1)
    img = ImageEnhance.Contrast(img).enhance(1)
    """plt.imshow(img)
    plt.show()"""

    shift = -110
    top, bottom = 340, 1200
    left, right = top + shift, bottom + shift

    im1 = img.crop((left, top, right, bottom))
    im1 = im1.resize((1600, 1600))

    im1.save(img_path.replace(".png", "_cropped2.png"))
    # im1.show()
