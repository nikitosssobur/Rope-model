from pydantic import BaseModel, field_validator, model_validator
from typing import Union



class ArgsValidator(BaseModel):
    x_coords: list[Union[float, int]]    
    y_coords: list[Union[float, int]] 
    radius: float
    moving_points_num: int
    moving_points_radius: float
    moving_points_velocities: list[Union[float, int]] 
    rope_mass: float
    tension_force: float
    stiffness: float
    viscosity: float
    rest_len: float
    point_num: int
    weighted_point_mass: float 
    weighted_point_radius: float
    ver_rope_point_radius: float
    ver_rope_mass: float
    ver_rope_stiffness: float
    ver_rope_viscosity: float
    attach_point_id: Union[int, None]
    surface_color: list[int]
    display_width: int
    display_height: int
    print_system_info: bool
    real_meter_dist: float
    fall_height: float
    sec_num: int
    time_steps_per_sec: int
    solver_method: str 
    rope_fragments_indexes: Union[list[int], None]


    @field_validator('x_coords')
    def check_x_coords(cls, value):
        if len(value) != 2:
            raise ValueError('x coords list must have only 2 values')
        if value[0] <= 0 or value[1] <= 0:
            raise ValueError('Fixed points x coords must be positive floats!')
        if value[0] >= value[1]:
            raise ValueError('Left fixed point x coord must be greater than right fixed point x coord!')
        return value
    
    
    @field_validator('y_coords')
    def check_y_coords(cls, value):
        if len(value) != 2:
            raise ValueError('y coords list must have only 2 values')
        if value[0] <= 0 or value[1] <= 0:
            raise ValueError('Fixed points y coords must be positive floats!')
        return value 
    

    @field_validator('radius', 'moving_points_radius', 'weighted_point_radius', 'ver_rope_point_radius')    
    def check_radius(cls, value):
        if value <= 0:
            raise ValueError('Fixed points, moving points, weighted point radius must be a positive float or integer!')
        return value 
    

    @field_validator('moving_points_num', 'point_num') 
    def check_moving_points_num(cls, value):
        if value <= 0:
            raise ValueError('Moving points number must be a positive integer!')
        return value
    

    @field_validator('moving_points_velocities')
    def check_moving_points_velocities(cls, value):
        if len(value) != 2:
            raise ValueError('velocities vector must have only 2 values')
        return value 
    

    @field_validator('tension_force')
    def check_tension_force(cls, value):
        if value <= 0:
            raise ValueError('Initial tension force of the horizontal rope must be a positive float!')
        return value


    @field_validator('rope_mass', 'ver_rope_mass', 'stiffness', 'viscosity', 'ver_rope_stiffness', 'ver_rope_viscosity')
    def check_rope_characteristic(cls, value):
        if value <= 0:
            raise ValueError('Ropes characteristics (mass, stiffness, viscosity) must be positive floats!')
        return value 
    

    @model_validator(mode="after")
    def check_attach_point_id(self) -> 'ArgsValidator':
        if isinstance(self.attach_point_id, Union[None]):
            self.attach_point_id = self.moving_points_num // 2 
        if not 0 <= self.attach_point_id <= self.moving_points_num - 1:
            raise ValueError('Attachment point index (attach_point_id) must be a positive integer in the segment: [0, horizontal rope points number]!')
        return self


    @field_validator('surface_color')
    def check_surface_color(cls, value):
        if len(value) != 3:
            raise ValueError('Color description list must have only 3 values (RGB format)')
        if not all(list(map(lambda v: 0 <= v <= 255, value))):
            raise ValueError('Color must be in RGB format: each value is integer and in segment [0, 255]!')
        return value        


    @field_validator('display_width', 'display_height')
    def check_display_data(cls, value):
        if value <= 0:
            raise ValueError('Display width and height must be positive integers!')
        return value


    @field_validator('real_meter_dist')
    def check_real_meter_dist(cls, value):
        if value <= 0:
            raise ValueError('Real meter distance must be a positive float!')
        return value


    @field_validator('fall_height')
    def check_fall_height(cls, value):
        if value <= 0:
            raise ValueError('Fall height must be a positive float!')
        return value


    @field_validator('sec_num')
    def check_sec_num(cls, value):
        if value <= 0:
            raise ValueError('Number of seconds (duration of the experiment) must be a positive integer!')
        return value 


    @field_validator('time_steps_per_sec')
    def check_time_steps_per_sec(cls, value):
        if value <= 0:
            raise ValueError('Number of time steps per one second (frames per second) must be a positive integer!')
        return value 


    @field_validator('solver_method')
    def check_solver_method(cls, value):
        methods = ['RK45', 'RK23','DOP853','Radau','BDF','LSODA']
        if value not in methods:
            raise ValueError(f'Inappropriate method value! Solver method value must be one of available strings: {methods}')
        return value
    

    @model_validator(mode="after")
    def check_rope_fragments_indexes(self) -> 'ArgsValidator':
        if isinstance(self.rope_fragments_indexes, list):
            all_ints = all([isinstance(val, int) for val in self.rope_fragments_indexes])
            if not all_ints:
                raise ValueError('List of rope fragments indexes must contain only integers!')
            
            bounds_check = all([-self.moving_points_num <= indx_val < self.moving_points_num for indx_val in self.rope_fragments_indexes])
            
            if not bounds_check:
                raise ValueError('All indexes in the rope fragments indexes list must be integers in segment [-hor_points_num, hor_points_num-1]')
       
        return self     
        