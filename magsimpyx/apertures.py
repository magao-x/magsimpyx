import numpy as np
from hcipy import *

def make_magaox_bump_mask(normalized=False, with_spiders=True):
    '''Make the Magellan bump mask.

    Parameters
    ----------
    normalized : boolean
        If this is True, the outer diameter will be scaled to 1. Otherwise, the
        diameter of the pupil will be 6.5 meters.
    with_spiders: boolean
        If this is False, the spiders will be left out.

    Returns
    -------
    Field generator
        The Magellan aperture.
    '''

    # TODO: Magnify bump mask vals. and preserve default vals

    magnification_factor = 6.5/9e-3 # Mag factor to scale 9 mm bump mask up to 6.5 m pupil diameter
    mask_inner = 2.79e-3 * magnification_factor # meter
    mask_outer = 8.604e-3 * magnification_factor # meter

    bump_mask_diameter = 0.5742e-3 * magnification_factor

    bump_mask_pos = [2.853e-3 * magnification_factor, -0.6705e-3 * magnification_factor] 

    radius = np.hypot(bump_mask_pos[0], bump_mask_pos[1])
    theta = np.arctan2(bump_mask_pos[1], bump_mask_pos[0]) - np.rad2deg(38.7747) + np.pi/2 # Adjusted bump angle to better center it on spider
    bump_mask_pos = [radius * np.cos(theta), radius * np.sin(theta)]

    pupil_diameter = 6.5 # meter
    spider_width1 = 0.1917e-3 * magnification_factor # meter 
    spider_width2 = 0.1917e-3  * magnification_factor # meter
    central_obscuration_ratio = mask_inner / mask_outer 
    spider_offset = [0, 0.34]  # meter

    if normalized:
        spider_width1 /= pupil_diameter
        spider_width2 /= pupil_diameter
        spider_offset = [x / pupil_diameter for x in spider_offset]
        bump_mask_pos = [x / pupil_diameter for x in bump_mask_pos]
        bump_mask_diameter /= pupil_diameter
        pupil_diameter = 1.0

    spider_offset = np.array(spider_offset)

    mirror_edge1 = (pupil_diameter / (2 * np.sqrt(2)), pupil_diameter / (2 * np.sqrt(2)))
    mirror_edge2 = (-pupil_diameter / (2 * np.sqrt(2)), pupil_diameter / (2 * np.sqrt(2)))
    mirror_edge3 = (pupil_diameter / (2 * np.sqrt(2)), -pupil_diameter / (2 * np.sqrt(2)))
    mirror_edge4 = (-pupil_diameter / (2 * np.sqrt(2)), -pupil_diameter / (2 * np.sqrt(2)))

    obstructed_aperture = make_obstructed_circular_aperture(mask_outer, central_obscuration_ratio)
    bump_mask = make_circular_aperture(bump_mask_diameter, center=bump_mask_pos) # Generate bump cover for Magellan pupil
    
    if not with_spiders:
        return obstructed_aperture

    spider1 = make_spider(spider_offset, mirror_edge1, spider_width1)
    spider2 = make_spider(spider_offset, mirror_edge2, spider_width1)
    spider3 = make_spider(-spider_offset, mirror_edge3, spider_width2)
    spider4 = make_spider(-spider_offset, mirror_edge4, spider_width2)

    def func(grid):
        return obstructed_aperture(grid) * spider1(grid) * spider2(grid) * spider3(grid) * spider4(grid) * (1 - bump_mask(grid))
    
    return func

def make_magaox_large_lyot_stop(normalized=False, with_spiders=True):
    '''Make the Magellan bump mask.

    Parameters
    ----------
    normalized : boolean
        If this is True, the outer diameter will be scaled to 1. Otherwise, the
        diameter of the pupil will be 6.5 meters.
    with_spiders: boolean
        If this is False, the spiders will be left out.

    Returns
    -------
    Field generator
        The Magellan aperture.
    '''

    # TODO: Magnify bump mask vals. and preserve default vals

    magnification_factor = 6.5/9e-3 # Mag factor to scale 9 mm bump mask up to 6.5 m pupil diameter
    mask_inner = 3.60017e-3 * magnification_factor # meter
    mask_outer = 8.02356e-3 * magnification_factor # meter

    bump_mask_diameter = 1.149e-3 * magnification_factor

    bump_mask_pos = [2.853e-3 * magnification_factor, -0.6705e-3 * magnification_factor] 

    radius = np.hypot(bump_mask_pos[0], bump_mask_pos[1])
    theta = np.arctan2(bump_mask_pos[1], bump_mask_pos[0]) - np.rad2deg(38.7747) + np.pi/2 # Adjusted bump angle to better center it on spider
    bump_mask_pos = [radius * np.cos(theta), radius * np.sin(theta)]

    pupil_diameter = 6.5 # meter
    spider_width1 = 0.3830e-3 * magnification_factor # meter 
    spider_width2 = 0.3830e-3  * magnification_factor # meter
    central_obscuration_ratio = mask_inner / mask_outer 
    spider_offset = [0, 0.34]  # meter

    if normalized:
        spider_width1 /= pupil_diameter
        spider_width2 /= pupil_diameter
        spider_offset = [x / pupil_diameter for x in spider_offset]
        bump_mask_pos = [x / pupil_diameter for x in bump_mask_pos]
        bump_mask_diameter /= pupil_diameter
        pupil_diameter = 1.0

    spider_offset = np.array(spider_offset)

    mirror_edge1 = (pupil_diameter / (2 * np.sqrt(2)), pupil_diameter / (2 * np.sqrt(2)))
    mirror_edge2 = (-pupil_diameter / (2 * np.sqrt(2)), pupil_diameter / (2 * np.sqrt(2)))
    mirror_edge3 = (pupil_diameter / (2 * np.sqrt(2)), -pupil_diameter / (2 * np.sqrt(2)))
    mirror_edge4 = (-pupil_diameter / (2 * np.sqrt(2)), -pupil_diameter / (2 * np.sqrt(2)))

    obstructed_aperture = make_obstructed_circular_aperture(mask_outer, central_obscuration_ratio)
    bump_mask = make_circular_aperture(bump_mask_diameter, center=bump_mask_pos) # Generate bump cover for Magellan pupil
    
    if not with_spiders:
        return obstructed_aperture

    spider1 = make_spider(spider_offset, mirror_edge1, spider_width1)
    spider2 = make_spider(spider_offset, mirror_edge2, spider_width1)
    spider3 = make_spider(-spider_offset, mirror_edge3, spider_width2)
    spider4 = make_spider(-spider_offset, mirror_edge4, spider_width2)

    def func(grid):
        return obstructed_aperture(grid) * spider1(grid) * spider2(grid) * spider3(grid) * spider4(grid) * (1 - bump_mask(grid))
    
    return func