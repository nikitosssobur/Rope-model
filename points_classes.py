import pygame as pg
import numpy as np
from typing import Union




class Point:
    '''
    Base point class. FixedPoint and MassPoint classes inherits methods from this class
    '''
    
    def __init__(self):
        self.point_color = 'black' 
        self.point_data = {'coordinates': [], 'radius': None}
    

    def draw_point(self):
        pg.draw.circle(self.surface, self.point_color, 
                       (self.x, self.y), self.radius)


    def set_surface(self, surface: pg.Surface):
        self.surface = surface


    def set_color(self, color):
        self.point_color = color    


    def set_radius(self, radius: Union[float, int]):
        self.radius = radius 
        self.point_data['radius'] = self.radius


    def set_coords(self, coords: Union[tuple, list, np.ndarray]):
        self.x, self.y = coords
        self.point_data['coordinates'] = np.array([self.x, self.y], dtype = 'd')



class MassPoint(Point):
    '''
        Class for weighted moving point of the horizontal and vertical ropes   
    '''
    
    def __init__(self, id: int, mass: Union[float, int], radius: Union[float, int], 
                 start_coords: Union[tuple, list, np.ndarray], 
                 start_velocity: Union[tuple, list, np.ndarray]):
        super().__init__()
        self.id = id
        self.set_coords(start_coords)
        self.vel_x, self.vel_y = start_velocity
        self.mass = mass
        self.set_radius(radius)
        self.point_color = 'black'
        self.point_data = {'id': self.id,
                        'coordinates': np.array([self.x, self.y]), 
                        'velocities': np.array([self.vel_x, self.vel_y]), 
                        'radius': self.radius,
                        'mass': self.mass}

        self.si_point_data = {'coordinates': None, 'velocities': None}


    def update_coords(self, updated_x: Union[float, int], updated_y: Union[float, int]):
        self.si_point_data['coordinates'] = np.array([updated_x, updated_y])
        self.x, self.y = self.meters_to_pixels()
        self.point_data['coordinates'] = np.array([self.x, self.y])

    
    def update_velocities(self, updated_vel_x: Union[float, int], updated_vel_y: Union[float, int]):
        self.si_point_data['velocities'] = np.array([updated_vel_x, updated_vel_y])
        self.vel_x, self.vel_y = self.meters_to_pixels(mode = 'vel')
        self.point_data['velocities'] = np.array([self.vel_x, self.vel_y])

    
    def get_data(self) -> dict: return self.point_data


    def set_mass(self, mass: Union[float, int]):
        if mass > 0:
            self.mass = mass
            self.point_data['mass'] = self.mass
        else:
            raise ValueError('Mass value of the point must be a positive real number!')
        

    def meters_to_pixels(self, mode: str = 'coord') -> np.ndarray:
        if mode == 'coord':
            pix_coords = self.si_point_data['coordinates'] * self.scale
            pix_coords += self.displace_vector
        elif mode == 'vel':
            pix_coords = self.si_point_data['velocities'] * self.scale
        return pix_coords
        
    
    def set_convert_data(self, scale: Union[float, int], displace_vector: np.ndarray):
        self.scale = scale
        self.displace_vector = displace_vector if isinstance(displace_vector, np.ndarray) else np.array(displace_vector)
        if list(self.si_point_data.values()) == [None, None]:
            self.si_point_data['coordinates'] = (self.point_data['coordinates'] - self.displace_vector) / self.scale
            self.si_point_data['velocities'] = self.point_data['velocities'] / self.scale       



class FixedPoint(Point):
    '''
        Class for fixed points on the edges! 
    '''
    def __init__(self, x_coord: Union[float, int], y_coord: Union[float, int], radius: Union[float, int], point_color = 'black'):
        self.radius = radius
        self.point_color = point_color
        self.__velocities = np.array([0, 0])
        self.point_data = {'coordinates': [],
                           'velocities': self.__velocities, 
                           'radius': self.radius,
                           'color': self.point_color}
        self.set_coords((x_coord, y_coord))
        self.si_point_data = {'coordinates': None, 'velocities': self.__velocities}


    def set_convert_data(self, scale: Union[float, int], displace_vector: np.ndarray):
        self.scale, self.displace_vector = scale, displace_vector
        if self.si_point_data['coordinates'] is None:
            self.si_point_data['coordinates'] = (self.point_data['coordinates'] - self.displace_vector) / self.scale