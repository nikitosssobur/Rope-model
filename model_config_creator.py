import json 
from rope_model import DEFAULT_FALL_HEIGHT, DEFAULT_REAL_METER_DIST


'''
THIS IS THE CODE TO CREATE AND SAVE THE DEFAULT MODEL CONFIG.
IT IS STRONGLY RECOMMENDED NOT TO CHANGE THIS JSON CONFIGURATION FILE! 
YOU CAN USE IT AS AN EXAMPLE OF MODEL SETUP.
'''



DEFAULT_CONFIG_PATH = r"./default_model_config.json" 


default_model_parameters = {'x_coords': [100, 750], 
                            'y_coords': [20, 20], 
                            'radius': 8, 
                            'moving_points_num': 11, 
                            'moving_points_radius': 3, 
                            'moving_points_velocities': [0, 0],
                            'rope_mass': 12, 
                            'tension_force': 800,
                            'stiffness': 3000,
                            'viscosity': 2000, 
                            'rest_len': 25, 
                            'point_num': 3, 
                            'weighted_point_mass': 90, 
                            'weighted_point_radius': 8,
                            'ver_rope_point_radius': 3,
                            'ver_rope_mass': 12,
                            'ver_rope_stiffness': 2000,
                            'ver_rope_viscosity': 2000, 
                            'attach_point_id': None,
                            'surface_color': [255, 255, 255],
                            'display_width': 850,
                            'display_height': 600,
                            'print_system_info': True,
                            'real_meter_dist': DEFAULT_REAL_METER_DIST, 
                            'fall_height': DEFAULT_FALL_HEIGHT, 
                            'sec_num': 35, 
                            'time_steps_per_sec': 500, 
                            'solver_method': 'Radau',
                            'rope_fragments_indexes': None}



with open(DEFAULT_CONFIG_PATH, "w") as f:
    json.dump(default_model_parameters, f, indent=4)
