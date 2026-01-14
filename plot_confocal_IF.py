import glob
import os
from readlif.reader import LifFile
from skimage import data
import napari
from tkinter import *
import cv2
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from PIL import ImageTk, Image, ImageEnhance
import time
from scipy.ndimage import rotate


class ZStackViewer:
    def __init__(self, root, images):
        self.root = root
        self.canvas = Canvas(root, width=600, height=600)
        self.canvas.pack()
        b = ImageTk.PhotoImage(images[0])
        self.canvas.create_image(400, 400, anchor=CENTER, image=b)

    def scroll_through_channels(self):
        pass


def create_rgb_old(gray, channel="red"):
    # print("create_rgb_old", end=": ")
    channel_index2 = 0
    if channel == "red":
        channel_index = 0
    elif channel == "green":
        channel_index = 1
    elif channel == "blue":
        channel_index = 2
    elif channel == "magenta":
        channel_index = 0
        channel_index2 = 2
    elif channel == "yellow":
        channel_index = 0
        channel_index2 = 1
    elif channel == "cyan":
        channel_index = 1
        channel_index2 = 2
    else:
        print("Please define right channel")
        return gray

    colored_image = list()
    for y in gray:
        b_temp = list()
        for x in y:
            temp_px = [0, 0, 0]
            temp_px[channel_index] = x
            if color_blind:
                temp_px[channel_index2] = x
            # print(temp_px)
            b_temp.append(temp_px)
        colored_image.append(b_temp)
    # colored_image = np.asarray(colored_image)
    return colored_image


def create_rgb(gray, channel="red"):
    # print("create_rgb", end=": ")
    """
    converts a gray image with shape (y, x, 1)

    :param gray:
    :param channel:
    :return: colorized_image: array-like
    """

    if channel == "red":
        return [[[0, 0, x] for x in y] for y in gray]
    elif channel == "green":
        return [[[0, x, 0] for x in y] for y in gray]
    elif channel == "blue":
        return [[[x, 0, 0] for x in y] for y in gray]
    elif channel == "magenta":
        return [[[x, 0, x] for x in y] for y in gray]
    elif channel == "yellow":
        return [[[x, x, 0] for x in y] for y in gray]
    elif channel == "cyan":
        return [[[0, x, x] for x in y] for y in gray]
    else:
        print("Please define right channel")
        return gray


def create_3d_plot():
    print(leica_project.get_image())
    image = leica_project.get_image(image_series)

    cmap_blue = colors.LinearSegmentedColormap.from_list("", ["black", "blue", "blue", "blue"])
    cmap_green = colors.LinearSegmentedColormap.from_list("", ["black", "green", "green", "green"])
    cmap_red = colors.LinearSegmentedColormap.from_list("", ["black", "red", "red", "red"])

    color_maps = [cmap_blue, cmap_green, cmap_red]
    channel_colors = ["blue", "green", "red"]
    if color_blind:
        channel_colors = ["cyan", "yellow", "magenta"]  # color-blind version
    thresholds = [35, 35, 65]
    alphas = [1, 1, 1]
    img_size = 256
    shift_factor = [0, img_size * 1.1, img_size * 2 * 1.1]

    fig = plt.figure(figsize=(16, 16))
    plt.subplots_adjust(left=0, right=1, bottom=0, top=1)
    ax = fig.gca(projection='3d')
    channel_rgbs = {0: [], 1: [], 2: []}

    for channel in range(0, 3):
        num_stacks = 0
        x = []
        y = []
        z = []
        c = []
        print("channel", channel)
        for idx, plane in enumerate(image.get_iter_z(c=channel)):
            num_stacks += 1
            # plane = image.get_plane(requested_dims={3: i})
            plane = plane.resize((img_size, img_size))
            plane = ImageEnhance.Contrast(plane).enhance(enhance_factor)
            # print(plane)
            # print(np.max(plane))
            a = np.asarray(plane)
            start = time.perf_counter()
            b = create_rgb(a, channel=channel_colors[channel])
            print(time.perf_counter() - start)

            start = time.perf_counter()
            b = create_rgb_old(a, channel=channel_colors[channel])
            print(time.perf_counter() - start)
            quit()
            channel_rgbs[channel].append(b)
            # print(b)
            """fig, ax = plt.subplots(1, 1)
            # ax.imshow(a, cmap=color_maps[channel])
            ax.imshow(b)
            plt.show()"""

            # kernel = np.ones((3, 3), np.uint8)
            # a = cv2.dilate(a, kernel, iterations=1)
            # a = cv2.erode(a, kernel, iterations=1)

            # print(a)
            roi = np.where(a > thresholds[channel])
            color_vals = a[roi].flatten()

            c.extend(color_vals)
            x.extend([x + shift_factor[channel] for x in roi[1]])
            y.extend([x + shift_factor[channel] for x in roi[0]])
            z.extend([idx] * len(color_vals))

        plt.style.use('dark_background')
        # plt.grid(color='k', linestyle='-', linewidth=0)
        ax.w_xaxis.pane.fill = False
        ax.w_yaxis.pane.fill = False
        ax.w_zaxis.pane.fill = False
        ax.set_facecolor('black')
        ax.grid(False)
        squeeze_factor = 1
        # ax.set_zlim(-num_stacks * squeeze_factor, num_stacks * squeeze_factor)
        print(len(x), len(y), len(z), len(c))
        # print(x, y, z, c)
        ax.scatter3D(x, y, z, c=c, cmap=color_maps[channel], s=5, alpha=alphas[channel], marker=",", edgecolors=None)

    plt.savefig(f"{image_series}/{image_series}_3d_plot.png")
    # plt.show()
    plt.close("all")


def create_z_stack_series_and():
    plt.style.use('default')
    print(leica_project.get_image())
    image = leica_project.get_image(image_series)
    print(image)
    channel_colors = ["blue", "green", "red", "white"]  # white = merged
    if color_blind:
        channel_colors = ["cyan", "yellow", "magenta", "white"]  # color-blind version
    channel_rgbs = {0: [], 1: [], 2: []}
    channel_superprojection = {0: [], 1: [], 2: []}

    for channel in range(0, 3):
        num_stacks = 0
        print("Process channel", channel)
        for idx, plane in enumerate(image.get_iter_z(c=channel)):
            if store_stacks_as_tif:
                plane.save(fr"{save_in}/{image_series}/{image_series}_CH{channel}_stack{idx}.tif")
            num_stacks += 1

            if channel == 1:
                plane = ImageEnhance.Brightness(plane).enhance(channel2_brightness)
                plane = ImageEnhance.Contrast(plane).enhance(channel2_contrast)
            if channel == 2:
                plane = ImageEnhance.Brightness(plane).enhance(channel3_brightness)
                plane = ImageEnhance.Contrast(plane).enhance(channel3_contrast)

            channel_superprojection[channel].append(np.asarray(plane))
            shift = 0
            top, bottom = 300, 900
            left, right = top + shift, bottom + shift

            # plane = plane.crop((left, top, right, bottom))
            # print(plane)
            # plane = plane.rotate(0, Image.Resampling.NEAREST, expand=0)
            a = np.asarray(plane)
            # start = time.perf_counter()
            b = create_rgb(a, channel=channel_colors[channel])
            # print(time.perf_counter() - start)
            channel_rgbs[channel].append(b)

    if super_projection:
        print("Create super-projections")
        for channel in range(0, 3):
            projection_temp = np.asarray(channel_superprojection[channel])
            projection = []
            for y in range(0, len(projection_temp[0])):
                projection.append(np.percentile(projection_temp[:, y], 96, axis=0))
            projection = np.asarray(projection).astype(int)
            colored_projection = create_rgb(projection, channel=channel_colors[channel])
            channel_rgbs[channel].append(colored_projection)
            # plt.imshow(colored_projection)
            # plt.show()

    # print(channel_rgbs)
    blues = channel_rgbs[0]
    greens = channel_rgbs[1]
    reds = channel_rgbs[2]

    print("Create merged image", end=" ")
    start = time.perf_counter()
    temp = np.asarray(list(zip(blues, greens, reds)))
    # print(temp[:, :, :, :, 0].shape)

    # rgb coloring
    px1 = np.max(temp[:, :, :, :, 0], axis=1)  # red channel
    px2 = np.max(temp[:, :, :, :, 1], axis=1)  # green channel
    px3 = np.max(temp[:, :, :, :, 2], axis=1)  # blue channel
    z_merged = [np.dstack((px1[stack], px2[stack], px3[stack])).reshape(1024, 1024, 3) for stack in range(0, len(px1))]
    print(time.perf_counter() - start, " sec")

    """start = time.perf_counter()
    z_merged = list()
    for b, g, r in zip(blues, greens, reds):
        merged = list()
        for y in range(len(b)):
            temp = list()
            for x in range(len(b[0])):
                if not color_blind:
                    temp.append([r[y][x][0], g[y][x][1], b[y][x][2]])
                else:
                    new_vals = [max([r[y][x][0], g[y][x][0], b[y][x][0]]),
                                max([r[y][x][1], g[y][x][1], b[y][x][1]]),
                                max([r[y][x][2], g[y][x][2], b[y][x][2]])]
                    temp.append(new_vals)
            merged.append(temp)
        z_merged.append(np.asarray(merged))
    print(time.perf_counter() - start)"""
    all_img = [blues, greens, reds, z_merged]

    print("Plot series of images")
    plt.figure(figsize=(len(z_merged) * 3, 12))
    plt.subplots_adjust(left=0.02, right=0.98, top=.98, bottom=0.02, wspace=0.02, hspace=0.02)
    col = 0
    row = 0
    channel_names_temp = list(channel_names.values())
    x_pos = int(len(blues[0][0]) * 0.05)
    y_pos = int(len(blues[0]) * 0.05)
    for idx, stack in enumerate(all_img):
        for idx2, img in enumerate(stack):

            if not 12 > idx2 > 5:
                pass
            plt.subplot2grid((4, len(z_merged)), (row, col))
            # cut image
            # img = img.rotate(18, Image.NEAREST, expand=0)
            # img.rotate(int(float(18)))
            img = np.asarray(img)
            if cut_image:
                img = img[100:950, 100:950]
            plt.imshow(img)
            if idx2 == 0:
                color = channel_colors[idx]
                plt.text(x=x_pos, y=y_pos, s=channel_names_temp[idx], color=color, va="top",
                         fontdict={"size": 20, "weight": "bold"})
            plt.axis("off")
            col += 1
        row += 1
        col = 0
    #plt.show()
    plt.savefig(f"{save_in}/{image_series}/{image_series}_z_stack_{image_series}_c_blind_{color_blind}_proj_{super_projection}.png")
    plt.savefig(f"{save_in}/{image_series}/{image_series}_z_stack_{image_series}_c_blind_{color_blind}_proj_{super_projection}.pdf")
    # plt.savefig(f"{image_series}/{image_series}_z_stack_{image_series}.tif")
    # plt.savefig(f"{image_series}/{image_series}_z_stack_{image_series}.pdf")
    # plt.show()
    plt.close("all")
    print("Plot save super-projection")
    for idx, images in enumerate(all_img):
        plt.figure(figsize=(16, 16), # dpi=600
                   )
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        plt.imshow(images[-1])
        plt.axis("off")
        channel_num = idx + 1
        plt.savefig(rf"{save_in}/{image_series}/super_projection_CH{channel_num}.png")
        plt.close()

    print("Save super-projection series")
    for idx, image in enumerate(z_merged[:-1]):
        plt.figure(figsize=(16, 16))
        plt.subplots_adjust(left=0, right=1, top=1, bottom=0)
        plt.imshow(image)
        plt.axis("off")
        stack_num = idx + 1
        plt.savefig(rf"{save_in}/{image_series}/merged_stack{stack_num}.png")
        plt.close()

    # print("create_gifs_from_images".replace("_", " "))
    # create_gifs_from_images()

    print("done")


def create_gifs_from_images():
    # create gifs and store in project subfolder
    img, *imgs = [Image.open(image) for image in sorted(glob.glob(rf"{save_in}/{image_series}/merged_stack*.png"),
                                                        key=os.path.getctime)]
    img.save(fp=rf"{save_in}/{image_series}/super_projection.gif",
             format='GIF', append_images=imgs,
             save_all=True,
             duration=int(950),
             # The display duration of each frame, in milliseconds.
             loop=10)


if __name__ == '__main__':

    # load Leica project
    project_name = "09032023_mSMC_CCNB1_Zc3_WT_WDH_Interphase.lif"
    leica_project = LifFile(rf"Raw\{project_name}")
    img_list = [i for i in leica_project.get_iter_image()]
    print(img_list)

    num_series = len(img_list)

    channel_names = {"blue": "DAPI",
                     "green": "Cyclin B1",
                     "red": "Tubulin",
                     "merge": "Merge"}

    rotation_angle = 0
    cut_image = False
    color_blind = True
    store_stacks_as_tif = False
    super_projection = True
    enhance_factor = 1.3

    save_in = project_name.replace(".lif", "")
    if not os.path.isdir(str(save_in)):
        os.mkdir(str(save_in))

    channel2_brightness = 1.2  # Green / Yellow
    channel2_contrast = 1.2   # max 1.5

    channel3_brightness = 1.2  # Red / Magenta
    channel3_contrast = 1.2  # max 1.5

    for i in range(0, num_series):
        if i != 10:
            pass

        print("\n", "Analyze Image: ", i)

        image_series = i
        save_name = f"{save_in}/{image_series}/{image_series}_z_stack_{image_series}_c_blind_{color_blind}_proj_{super_projection}.png"
        if os.path.isfile(save_name):
            print("Already analyzed")
            continue
        save_in = project_name.replace(".lif", "")
        if not os.path.isdir(str(save_in) + "/" + str(image_series)):
            os.mkdir(str(save_in) + "/" + str(image_series))
        create_z_stack_series_and()

        # create_gifs_from_images()

        # create_3d_plot()
    quit()
    root = Tk()
    app = ZStackViewer
    root.mainloop()
