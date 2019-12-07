import random

import numpy as np
from genevis.render import RaycastRenderer
from genevis.transfer_function import TFColor
from volume.volume import GradientVolume, Volume
from collections.abc import ValuesView
from scipy.interpolate import RegularGridInterpolator
from numpy import linspace, zeros, array
import math




def TriLinear(volume: Volume, x: float, y: float, z: float):
    x0 = int(math.floor(x))
    x1 = math.ceil(x)
    y0 = int(math.floor(y))
    y1 = math.ceil(y)
    z0 = int(math.floor(z))
    z1 = math.ceil(z)

    if x1 < 0 or y1 < 0 or z1 < 0 or x1 >= volume.dim_x or y1 >= volume.dim_y or z1 >= volume.dim_z:
        return 0

    alpha = (x - x0)
    beta = (y - y0)
    teta = (z - z0)

    c000 = volume.data[x0, y0, z0]
    c100 = volume.data[x1, y0, z0]
    c001 = volume.data[x0, y0, z1]
    c101 = volume.data[x1, y0, z1]
    c010 = volume.data[x0, y1, z0]
    c110 = volume.data[x1, y1, z0]
    c011 = volume.data[x0, y1, z1]
    c111 = volume.data[x1, y1, z1]

    c00=c000*(1-alpha)+c100*(alpha)
    c01=c001*(1-alpha)+c101*(alpha)
    c10=c010*(1-alpha)+c110*(alpha)
    c11=c011*(1-alpha)+c111*(alpha)

    c0 = c00*(1 - beta) + c10*(beta)
    c1 = c01*(1 - beta) + c11*(beta)

    point = c0*(1 - teta) + c1*(teta)
    return point
# TODO: Implement trilinear interpolation
def get_voxel(volume: Volume, x: float, y: float, z: float):
    """
    Retrieves the value of a voxel for the given coordinates.
    :param volume: Volume from which the voxel will be retrieved.
    :param x: X coordinate of the voxel
    :param y: Y coordinate of the voxel
    :param z: Z coordinate of the voxel
    :return: Voxel value
    """
    if x < 0 or y < 0 or z < 0 or x >= volume.dim_x or y >= volume.dim_y or z >= volume.dim_z:
        return 0

    x = int(math.floor(x))
    y = int(math.floor(y))
    z = int(math.floor(z))

    return volume.data[x, y, z]

def mip_get_voxel(volume: Volume, x: float, y: float, z: float):
    """
    Retrieves the value of a voxel for the given coordinates.
    :param volume: Volume from which the voxel will be retrieved.
    :param x: X coordinate of the voxel
    :param y: Y coordinate of the voxel
    :param z: Z coordinate of the voxel
    :return: Voxel value
    """
    if x < 0 or y < 0 or z < 0 or x >= volume.dim_x or y >= volume.dim_y or z >= volume.dim_z:
        return 0

    return TriLinear(volume,x,y,z)




class RaycastRendererImplementation(RaycastRenderer):
    """
    Class to be implemented.
    """

    def clear_image(self):
        """Clears the image data"""
        self.image.fill(0)

    # TODO: Implement trilinear interpolation


    def render_slicer(self, view_matrix: np.ndarray, volume: Volume, image_size: int, image: np.ndarray):
        # Clear the image
        self.clear_image()

        # U vector. See documentation in parent's class
        u_vector = view_matrix[0:3]
        print("u_vector",u_vector)
        # V vector. See documentation in parent's class
        v_vector = view_matrix[4:7]
        print("v_vector",v_vector)
        # View vector. See documentation in parent's class
        view_vector = view_matrix[8:11]
        print ("view_vector",view_vector)
        # Center of the image. Image is squared
        image_center = image_size / 2
        print ("image_size",image_size,"image_center",image_center)
        # Center of the volume (3-dimensional)
        volume_center = [volume.dim_x / 2, volume.dim_y / 2, volume.dim_z / 2]
        volume_maximum = volume.get_maximum()
        print("volume center",volume_center,", volume_maximum",volume_maximum)
        # Define a step size to make the loop faster
        step = 2 if self.interactive_mode else 1

        for i in range(0, image_size, step):
            for j in range(0, image_size, step):
                # Get the voxel coordinate X
                voxel_coordinate_x = u_vector[0] * (i - image_center) + v_vector[0] * (j - image_center) + \
                                     volume_center[0]


                # Get the voxel coordinate Y
                voxel_coordinate_y = u_vector[1] * (i - image_center) + v_vector[1] * (j - image_center) + \
                                     volume_center[1]

                # Get the voxel coordinate Z
                voxel_coordinate_z = u_vector[2] * (i - image_center) + v_vector[2] * (j - image_center) + \
                                     volume_center[2]

                # Get voxel value
                value = get_voxel(volume, voxel_coordinate_x, voxel_coordinate_y, voxel_coordinate_z)

                # Normalize value to be between 0 and 1
                red = value / volume_maximum
                green = red
                blue = red
                alpha = 1.0 if red > 0 else 0.0

                # Compute the color value (0...255)
                red = math.floor(red * 255) if red < 255 else 255
                green = math.floor(green * 255) if green < 255 else 255
                blue = math.floor(blue * 255) if blue < 255 else 255
                alpha = math.floor(alpha * 255) if alpha < 255 else 255

                # Assign color to the pixel i, j
                image[(j * image_size + i) * 4] = red
                image[(j * image_size + i) * 4 + 1] = green
                image[(j * image_size + i) * 4 + 2] = blue
                image[(j * image_size + i) * 4 + 3] = alpha


    # TODO: Implement MIP function
    def render_mip(self , view_matrix: np.ndarray, volume: Volume, image_size: int, image: np.ndarray):
        self.clear_image()
        u_vector = view_matrix[0:3]
        v_vector = view_matrix[4:7]
        view_vector = view_matrix[8:11]
        image_center = image_size / 2
        diagonal=int(image_center)
        img_neg=(-1)*image_center
        volume_center = [volume.dim_x / 2, volume.dim_y / 2, volume.dim_z / 2]
        volume_maximum = volume.get_maximum()
        # step = 2 if self.interactive_mode else 1
        for i in range(0, image_size, 1):
            for j in range(0, image_size, 1):
                m = 0
                for k in range(int(img_neg), int(image_center) ,10):

                    voxel_x = u_vector[0] * (i - image_center) + v_vector[0] * (j - image_center) + \
                                         volume_center[0]+view_vector[0]*k
                    voxel_y = u_vector[1] * (i - image_center) + v_vector[1] * (j - image_center) + \
                                         volume_center[1] + view_vector[1] * k
                    voxel_z = u_vector[2] * (i - image_center) + v_vector[2] * (j - image_center) + \
                                         volume_center[2] + view_vector[2] * k
                    tot_vox = mip_get_voxel(volume,voxel_x,voxel_y,voxel_z)
                    if(tot_vox>m):
                        m=tot_vox
                red = m / volume_maximum
                green = red
                blue = red
                alpha = 1.0 if red > 0 else 0.0

                                # Compute the color value (0...255)
                red = math.floor(red * 255) if red < 255 else 255
                green = math.floor(green * 255) if green < 255 else 255
                blue = math.floor(blue * 255) if blue < 255 else 255
                alpha = math.floor(alpha * 255) if alpha < 255 else 255

                                # Assign color to the pixel i, j
                image[(j * image_size + i) * 4] = red
                image[(j * image_size + i) * 4 + 1] = green
                image[(j * image_size + i) * 4 + 2] = blue
                image[(j * image_size + i) * 4 + 3] = alpha


    # TODO: Implement Compositing function. TFColor is already imported. self.tfunc is the current transfer function.
    def render_compositing(self, view_matrix: np.ndarray, volume: Volume, image_size: int, image: np.ndarray):
        self.clear_image()
        u_vector = view_matrix[0:3]
        v_vector = view_matrix[4:7]
        view_vector = view_matrix[8:11]
        image_center = image_size / 2
        diagonal = int(image_center)
        neg_diagonal = -1 * diagonal
        img_neg = (-1) * image_center
        volume_center = [volume.dim_x / 2, volume.dim_y / 2, volume.dim_z / 2]
        volume_maximum = volume.get_maximum()
        # step = 2 if self.interactive_mode else 1
        self.tfunc.set_test_function()
        for i in range(0, image_size, 1):
            for j in range(0, image_size, 1):
                c_color=TFColor(0,0,0,1)
                for k in range(int(img_neg), int(image_center), 10):

                    voxel_x = u_vector[0] * (i - image_center) + v_vector[0] * (j - image_center) + \
                              volume_center[0] + view_vector[0] * k
                    voxel_y = u_vector[1] * (i - image_center) + v_vector[1] * (j - image_center) + \
                              volume_center[1] + view_vector[1] * k
                    voxel_z = u_vector[2] * (i - image_center) + v_vector[2] * (j - image_center) + \
                              volume_center[2] + view_vector[2] * k
                    # if (voxel_x < 0 or voxel_y >= volume.dim_x
                    #         or voxel_y  < 0 or voxel_y >= volume.dim_y
                    #         or voxel_z < 0 or voxel_z >= volume.dim_z):
                    #         continue

                    tot_vox = get_voxel(volume, voxel_x, voxel_y, voxel_z)
                    color=self.tfunc.get_color(math.floor(tot_vox))
                    c_color.r = ((color.a) * color.r) + ((1 - color.a) * c_color.r)
                    c_color.g = ((color.a) * color.g) + ((1 - color.a) * c_color.g)
                    c_color.b = ((color.a) * color.b) + ((1 - color.a) * c_color.b)


                red = c_color.r * 255 if c_color.r < 255 else 255
                green = c_color.g * 255 if c_color.g < 255 else 255
                blue =  c_color.b * 255 if  c_color.b < 255 else 255
                alpha = c_color.a * 255 if c_color.a < 255 else 255

                # Assign color to the pixel i, j
                image[(j * image_size + i) * 4] = red
                image[(j * image_size + i) * 4 + 1] = green
                image[(j * image_size + i) * 4 + 2] = blue
                image[(j * image_size + i) * 4 + 3] = alpha




    # TODO: Implement function to render multiple energy volumes and annotation volume as a silhouette.
    def render_mouse_brain(self, view_matrix: np.ndarray, annotation_volume: Volume, energy_volumes: dict,
                           image_size: int, image: np.ndarray):
        # TODO: Implement your code considering these volumes (annotation_volume, and energy_volumes)
        self.clear_image()
        # an_volume=GradientVolume(annotation_volume)
        # an_volume.compute()
        mouse_brain_annotation(self, view_matrix,annotation_volume,energy_volumes,image_size,image)
        # gradient_volume=GradientVolume(annotation_volume)
        # gradient_volume.compute()
        # colors=TFColor(120,120,120,0)
        # TF_render(self, view_matrix, image_size, image,gradient_volume,colors, annotation_volume)
        #
        # for key in energy_volumes:
        #     energy_colors=TFColor(random.random(),random.random(),random.random(),0)
        #     TF_render(self, view_matrix, image_size, image, energy_volumes[key],energy_colors, annotation_volume )


def mouse_brain_energy(self, an_x: float, an_y: float, an_z:float, view_matrix: np.ndarray, energy_volumes:dict, image_size: int, image: np.ndarray):
    u_vector = view_matrix[0:3]
    v_vector = view_matrix[4:7]
    view_vector = view_matrix[8:11]
    image_center = image_size / 2
    en=list(energy_volumes.values())[0]
    en_vol_center=[en.dim_x/2, en.dim_y/2, en.dim_z/2]
    en_diag=math.floor(math.sqrt(
            en.dim_x * en.dim_x\
            + en.dim_y * en.dim_y \
            + en.dim_z * en.dim_z))
    neg_en_diag=(-1)*en_diag
    threshold=15.0
    step=50
    for key in energy_volumes.values():
        for i in range(0, image_size, step):
            for j in range(0, image_size, step):

                c_color = TFColor(0, 0, 0, 1)
                for k in range(int(neg_en_diag), int(en_diag), step):
                    voxel_x = u_vector[0] * (i - image_center) + v_vector[0] * (j - image_center) + \
                              en_vol_center[0] + view_vector[0] * k
                    voxel_y = u_vector[1] * (i - image_center) + v_vector[1] * (j - image_center) + \
                              en_vol_center[1] + view_vector[1] * k
                    voxel_z = u_vector[2] * (i - image_center) + v_vector[2] * (j - image_center) + \
                              en_vol_center[2] + view_vector[2] * k
                    tot_vox = mip_get_voxel(key, voxel_x, voxel_y, voxel_z)

                    if((an_x>voxel_x and voxel_x>=0)
                            or (an_y>voxel_y and voxel_y>=0)
                                or (an_z>voxel_z and voxel_z>=0)) and tot_vox>threshold:
                        colors = TFColor(random.random(),random.random(),random.random(),random.random())
                        c_color.r = ((colors.a) * colors.r) + ((1 - colors.a) * c_color.r)
                        c_color.g = ((colors.a) * colors.g) + ((1 - colors.a) * c_color.g)
                        c_color.b = ((colors.a) * colors.b) + ((1 - colors.a) * c_color.b)
                        c_color.a = ((colors.a) * colors.a) + ((1 - colors.a) * c_color.a)

                red = math.floor(c_color.r * 255) if c_color.r < 255 else 255
                green = math.floor(c_color.g * 255) if c_color.g < 255 else 255
                blue = math.floor(c_color.b * 255) if c_color.b < 255 else 255
                alpha = math.floor(c_color.a * 255) if c_color.a < 255 else 255

                # Assign color to the pixel i, j
                image[(j * image_size + i) * 4] = red
                image[(j * image_size + i) * 4 + 1] = green
                image[(j * image_size + i) * 4 + 2] = blue
                image[(j * image_size + i) * 4 + 3] = alpha


def mouse_brain_annotation(self, view_matrix: np.ndarray, annotation_volume: Volume, energy_volumes:dict, image_size: int, image: np.ndarray):
    u_vector = view_matrix[0:3]
    v_vector = view_matrix[4:7]
    view_vector = view_matrix[8:11]
    image_center = image_size / 2
    an_vol_center= [annotation_volume.dim_x / 2, annotation_volume.dim_y / 2, annotation_volume.dim_z / 2]
    max_grad_mag= self.annotation_gradient_volume.get_max_gradient_magnitude()
    an_diag=math.floor(math.sqrt(
            annotation_volume.dim_x * annotation_volume.dim_x\
            + annotation_volume.dim_y * annotation_volume.dim_y \
            + annotation_volume.dim_z * annotation_volume.dim_z))

    neg_an_diag=(-1)*an_diag
    step = 100
    for i in range(0, image_size, step):
        for j in range(0, image_size, step):
            c_color = TFColor(0, 0, 0, 0)

            for k in range(int(neg_an_diag), int(an_diag), step):
                voxel_x = u_vector[0] * (i - image_center) + v_vector[0] * (j - image_center) + \
                          an_vol_center[0] + view_vector[0] * k
                voxel_y = u_vector[1] * (i - image_center) + v_vector[1] * (j - image_center) + \
                          an_vol_center[1] + view_vector[1] * k
                voxel_z = u_vector[2] * (i - image_center) + v_vector[2] * (j - image_center) + \
                          an_vol_center[2] + view_vector[2] * k
                tot_vox = mip_get_voxel(annotation_volume, voxel_x, voxel_y, voxel_z)
                if tot_vox>0:
                    grad_mag = self.annotation_gradient_volume.get_gradient(int(math.floor(voxel_x)),
                                                                                int(math.floor(voxel_y)),
                                                                                int(math.floor(voxel_z))).magnitude
                    grad=grad_mag / (max_grad_mag*1.75)
                    colors = TFColor(0.7, 0.4, 0.5, grad)
                    c_color.r = ((colors.a) * colors.r) + ((1 - colors.a) * c_color.r)
                    c_color.g = ((colors.a) * colors.g) + ((1 - colors.a) * c_color.g)
                    c_color.b = ((colors.a) * colors.b) + ((1 - colors.a) * c_color.b)
                    c_color.a = ((colors.a) * colors.a) + ((1 - colors.a) * c_color.a)

            red = math.floor(c_color.r * 255) if c_color.r < 255 else 255
            green = math.floor(c_color.g * 255) if c_color.g < 255 else 255
            blue = math.floor(c_color.b * 255) if c_color.b < 255 else 255
            alpha = math.floor(c_color.a * 255) if c_color.a < 255 else 255

            # Assign color to the pixel i, j
            image[(j * image_size + i) * 4] = red
            image[(j * image_size + i) * 4 + 1] = green
            image[(j * image_size + i) * 4 + 2] = blue
            image[(j * image_size + i) * 4 + 3] = alpha


    if len(energy_volumes) == 0:
        return

    mouse_brain_energy(self, annotation_volume.dim_x, annotation_volume.dim_y, annotation_volume.dim_z, view_matrix,
                           energy_volumes, image_size, image)



# def TF_render(self, view_matrix: np.ndarray, image_size: int,image: np.ndarray, gradient_volume: Volume,colors:TFColor):
#         u_vector = view_matrix[0:3]
#         v_vector = view_matrix[4:7]
#         view_vector = view_matrix[8:11]
#         image_center = image_size / 2
#         img_neg = (-1) * image_center
#         volume_center = [volume.dim_x / 2, volume.dim_y / 2, volume.dim_z / 2]
#         fv = random.random()
#         radius = 1
#         alpha_v = 1
#         step = 50
#         for i in range(0, image_size, step):
#             for j in range(0, image_size, step):
#                 c_color = TFColor(0, 0, 0, 1)
#                 for k in range(int(img_neg), int(image_center), step):
#                     voxel_x = u_vector[0] * (i - image_center) + v_vector[0] * (j - image_center) + \
#                               volume_center[0] + view_vector[0] * k
#                     voxel_y = u_vector[1] * (i - image_center) + v_vector[1] * (j - image_center) + \
#                               volume_center[1] + view_vector[1] * k
#                     voxel_z = u_vector[2] * (i - image_center) + v_vector[2] * (j - image_center) + \
#                               volume_center[2] + view_vector[2] * k
#                     fx = get_voxel(gradient_volume, voxel_x, voxel_y, voxel_z)
#
#                     if fx>0:
#
#                         var_gradient = self.annotation_gradient_volume.get_gradient( int(math.floor(voxel_x)), int(math.floor(voxel_y)), int(math.floor(voxel_z)).magnitude   )
#
#
#
#                     if mag_grad == 0 and fx == fv:
#                         alpha = alpha_v * 1
#
#                     elif mag_grad > 0 and fx - (radius * mag_grad) <= fv and fx + (radius * mag_grad) >= fv:
#                         alpha = alpha_v * (1 - 1 / radius * abs((fv - fx) / mag_grad))
#
#                     else:
#                         alpha = 0
#
#                     CTotal=CTotal*(1-alpha)+red*alpha
#                     Alpha=Alpha+(1-Alpha)*alpha
#
#
#                 red = c_color.r * 255 if c_color.r < 255 else 255
#                 green = c_color.g * 255 if c_color.g < 255 else 255
#                 blue = c_color.b * 255 if c_color.b < 255 else 255
#                 Alpha = c_color.a * 255 if c_color.a < 255 else 255
#
#         #Assign color to the pixel i, j
#                 image[(j * image_size + i) * 4] = red
#                 image[(j * image_size + i) * 4 + 1] = green
#                 image[(j * image_size + i) * 4 + 2] = blue
#                 image[(j * image_size + i) * 4 + 3] = Alpha
#
#



# class GradientVolumeImpl(GradientVolume):
#     # TODO: Implement gradient compute function. See parent class to check available attributes.
#     def compute(self):
#         pass
