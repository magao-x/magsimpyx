import numpy as np
from hcipy import make_obstructed_circular_aperture, make_circular_aperture, make_spider_infinite, Field, make_regular_polygon_aperture, make_elliptical_aperture

__all__ = [
    'make_magaox_bump_mask',
    'make_magaox_large_lyot_stop',
	'make_gmt_lyot_aperture'
]

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

    magnification_factor = 6.5/9e-3 # Mag factor to scale 9 mm bump mask up to 6.5 m pupil diameter
    mask_inner = 2.79e-3 * magnification_factor # meter
    mask_outer = 8.604e-3 * magnification_factor # meter

    bump_mask_diameter = 0.5742e-3 * magnification_factor

    bump_mask_pos = np.array([2.853e-3, -0.6705e-3]) * magnification_factor

    radius = np.hypot(bump_mask_pos[0], bump_mask_pos[1])
    theta = np.arctan2(bump_mask_pos[1], bump_mask_pos[0]) - np.deg2rad(25.2) #+ np.pi/2 # Adjusted bump angle to better center it on spider
    bump_mask_pos = [radius * np.cos(theta), radius * np.sin(theta)]

    pupil_diameter = 6.5  # meter
    spider_width1 = 0.1917e-3 * magnification_factor  # meter 
    spider_width2 = 0.1917e-3  * magnification_factor  # meter
    central_obscuration_ratio = mask_inner / mask_outer 
    spider_offset = np.array([0.0, 0.34])  # meter

    if normalized:
        spider_width1 /= pupil_diameter
        spider_width2 /= pupil_diameter
        spider_offset /= pupil_diameter
        bump_mask_pos /= pupil_diameter
        bump_mask_diameter /= pupil_diameter
        pupil_diameter = 1.0

    obstructed_aperture = make_obstructed_circular_aperture(mask_outer, central_obscuration_ratio)
    bump_mask = make_circular_aperture(bump_mask_diameter, center=bump_mask_pos)  # Generate bump cover for the MagAO-X DM
    
    if not with_spiders:
        return obstructed_aperture
	
	# spider offsets corrections based on bumpMask fits file from J. Males.
    spider1 = make_spider_infinite(spider_offset - np.array([0, -0.0185]), 45.0, spider_width1)
    spider2 = make_spider_infinite(-spider_offset - np.array([0, 0.020]), -45.0, spider_width1)
    spider3 = make_spider_infinite(-spider_offset + np.array([0, -0.022]), 45.0 + 180.0, spider_width2)
    spider4 = make_spider_infinite(spider_offset + np.array([0, 0.023]), -45.0 + 180.0, spider_width2)

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
    magnification_factor = 6.5/9e-3 # Mag factor to scale 9 mm bump mask up to 6.5 m pupil diameter
    mask_inner = 3.60017e-3 * magnification_factor # meter
    mask_outer = 8.02356e-3 * magnification_factor # meter

    bump_mask_diameter = 1.149e-3 * magnification_factor

    bump_mask_pos = np.array([2.853e-3, -0.6705e-3]) * magnification_factor

    radius = np.hypot(bump_mask_pos[0], bump_mask_pos[1])
    theta = np.arctan2(bump_mask_pos[1], bump_mask_pos[0]) - np.deg2rad(25.2)
    bump_mask_pos = [radius * np.cos(theta), radius * np.sin(theta)]

    pupil_diameter = 6.5 # meter
    spider_width1 = 0.3830e-3 * magnification_factor # meter 
    spider_width2 = 0.3830e-3  * magnification_factor # meter
    central_obscuration_ratio = mask_inner / mask_outer 
    spider_offset = np.array([0, 0.34])  # meter

    if normalized:
        spider_width1 /= pupil_diameter
        spider_width2 /= pupil_diameter
        spider_offset /= pupil_diameter
        bump_mask_pos /= pupil_diameter
        bump_mask_diameter /= pupil_diameter
        pupil_diameter = 1.0

    obstructed_aperture = make_obstructed_circular_aperture(mask_outer, central_obscuration_ratio)
    bump_mask = make_circular_aperture(bump_mask_diameter, center=bump_mask_pos) # Generate bump cover for Magellan pupil
    
    if not with_spiders:
        return obstructed_aperture

	# spider offsets corrections based on bumpMask fits file from J. Males.
    spider1 = make_spider_infinite(spider_offset - np.array([0, -0.0185]), 45.0, spider_width1)
    spider2 = make_spider_infinite(-spider_offset - np.array([0, 0.020]), -45.0, spider_width1)
    spider3 = make_spider_infinite(-spider_offset + np.array([0, -0.022]), 45.0 + 180.0, spider_width2)
    spider4 = make_spider_infinite(spider_offset + np.array([0, 0.023]), -45.0 + 180.0, spider_width2)

    def func(grid):
        return obstructed_aperture(grid) * spider1(grid) * spider2(grid) * spider3(grid) * spider4(grid) * (1 - bump_mask(grid))
    
    return func

def make_gmt_lyot_aperture(normalized=False, with_spiders=True, return_segments=False, undersize=1, spider_width=1, truss_oversize=1.0):
	'''Make the Giant Magellan Telescope aperture.

	The primary mirror parameters come from the GMT Observatory Architecture Documentation (GMT-REQ-03215, Rev. F):
	https://www.gmto.org/resources/slpdr/ . Small corrections have been applied to match to the actual pupil from internal GMTO files.

	Parameters
	----------
	normalized : boolean
		If this is True, the outer diameter will be scaled to 1. Otherwise, the
		diameter of the pupil will be 25.448 meters.
	with_spiders : boolean
		If this is False, the spiders will be left out. Default: True.
	return_segments : boolean
		If this is True, the segments will also be returned as a list of Field generators.

	Returns
	-------
	Field generator
		The GMT aperture.
	elt_segments : list of Field generators
		The segments. Only returned when `return_segments` is True.
	'''
	gmt_outer_diameter = 25.448
	segment_size = (8.365 - 0.072) * undersize
	off_axis_segment_size = (8.365 - 0.015) * undersize

	# The spider truss to hold the secondary
	spider_1_width = 0.119 * spider_width
	spider_2_width = 0.115 * spider_width
	radius_spider_1 = 2.386
	radius_spider_2 = 2.409
	offset_spider_2 = -0.05
	truss_size = 4.93 * truss_oversize

	# 0.359 mm from the off-axis segment to the on-axis segment
	segment_gap = 0.359 + 0.088
	off_axis_tilt = np.deg2rad(13.522)
	central_hole_size = 3.495 / segment_size
	segment_distance = (segment_size / 2 + off_axis_segment_size / 2 * np.cos(off_axis_tilt) + segment_gap) / undersize

	if normalized:
		segment_size /= gmt_outer_diameter
		off_axis_segment_size /= gmt_outer_diameter
		segment_gap /= gmt_outer_diameter
		segment_distance /= gmt_outer_diameter
		truss_size /= gmt_outer_diameter
		spider_1_width /= gmt_outer_diameter
		spider_2_width /= gmt_outer_diameter
		radius_spider_1 /= gmt_outer_diameter
		radius_spider_2 /= gmt_outer_diameter
		offset_spider_2 /= gmt_outer_diameter

	def make_diverging_spider(position, start_width, divergence, orientation):
		def func(grid):
			y = grid.shifted(position).rotated(orientation).y
			x = grid.shifted(position).rotated(orientation).x
			return Field(abs(y) < (np.sin(divergence) * abs(x) + start_width / 2) * (x >= 0), grid)
		return func

	def make_central_gmt_segment(grid):
		center_segment = make_obstructed_circular_aperture(segment_size, central_hole_size)(grid)
		if with_spiders:
			spider_attachement_mask = 1 - make_regular_polygon_aperture(3, truss_size, angle=np.pi, center=None)(grid)

			spider_mask = grid.ones()
			for i in range(3):
				offset_angle = 2 * np.pi / 3 * i
				spider_1_start = radius_spider_1 * np.array([np.sin(offset_angle), np.cos(offset_angle)])
				spider_mask *= 1 - make_diverging_spider(spider_1_start, spider_1_width, np.deg2rad(1.16), np.deg2rad(9.83) + offset_angle)(grid)
				spider_mask *= 1 - make_diverging_spider(spider_1_start, spider_1_width, np.deg2rad(1.16), np.deg2rad(180.0 - 9.83) + offset_angle)(grid)

				spider_2_start = radius_spider_2 * np.array([np.sin(offset_angle), np.cos(offset_angle)]) + offset_spider_2 * np.array([np.cos(offset_angle), -np.sin(offset_angle)])
				spider_mask *= 1 - make_diverging_spider(spider_2_start, spider_2_width, np.deg2rad(0.0), np.deg2rad(-11.2) + offset_angle)(grid)

			return center_segment * spider_mask * spider_attachement_mask
		else:
			return center_segment

	segment_functions = [make_central_gmt_segment]
	for i in range(6):
		rotation_angle = np.pi / 3 * i
		xc = segment_distance * np.cos(rotation_angle)
		yc = segment_distance * np.sin(rotation_angle)

		aperture = make_elliptical_aperture([off_axis_segment_size * np.cos(off_axis_tilt), off_axis_segment_size], center=[xc, yc], angle=-rotation_angle)
		segment_functions.append(aperture)

	# Center segment obscurations
	def make_aperture(grid):
		aperture = grid.zeros()
		for segment in segment_functions:
			aperture += segment(grid)
		return aperture

	if return_segments:
		return make_aperture, segment_functions
	else:
		return make_aperture