import pygame as pg
import numpy as np





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


    def set_surface(self, surface):
        self.surface = surface


    def set_color(self, color):
        self.point_color = color    


    def set_radius(self, radius):
        self.radius = radius 
        self.point_data['radius'] = self.radius


    def set_coords(self, coords):
        self.x, self.y = coords
        self.point_data['coordinates'] = np.array([self.x, self.y], dtype = 'd')



class MassPoint(Point):
    '''
        Class for weighted moving point of the horizontal and vertical ropes   
    '''
    
    def __init__(self, id, mass, radius, start_coords, start_velocity):
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


    def update_coords(self, updated_x, updated_y):
        self.si_point_data['coordinates'] = np.array([updated_x, updated_y])
        self.x, self.y = self.meters_to_pixels()
        self.point_data['coordinates'] = np.array([self.x, self.y])

    
    def update_velocities(self, updated_vel_x, updated_vel_y):
        self.si_point_data['velocities'] = np.array([updated_vel_x, updated_vel_y])
        self.vel_x, self.vel_y = self.meters_to_pixels(mode = 'vel')
        self.point_data['velocities'] = np.array([self.vel_x, self.vel_y])

    
    def get_data(self): return self.point_data


    def set_mass(self, mass):
        if mass > 0:
            self.mass = mass
            self.point_data['mass'] = self.mass
        else:
            raise ValueError('Mass value of the point must be a positive real number!')
        

    def meters_to_pixels(self, mode = 'coord'):
        if mode == 'coord':
            pix_coords = self.si_point_data['coordinates'] * self.scale
            pix_coords += self.displace_vector
        elif mode == 'vel':
            pix_coords = self.si_point_data['velocities'] * self.scale
        return pix_coords
        
    
    def set_convert_data(self, scale, displace_vector):
        self.scale = scale
        self.displace_vector = displace_vector
        if list(self.si_point_data.values()) == [None, None]:
            self.si_point_data['coordinates'] = (self.point_data['coordinates'] - self.displace_vector) / self.scale
            self.si_point_data['velocities'] = self.point_data['velocities'] / self.scale       



class FixedPoint(Point):
    '''
        Class for fixed points on the edges! 
    '''
    def __init__(self, x_coord, y_coord, radius, point_color = 'black'):
        self.radius = radius
        self.point_color = point_color
        self.__velocities = np.array([0, 0])
        self.point_data = {'coordinates': [],
                           'velocities': self.__velocities, 
                           'radius': self.radius,
                           'color': self.point_color}
        self.set_coords((x_coord, y_coord))
        self.si_point_data = {'coordinates': None, 'velocities': self.__velocities}


    def set_convert_data(self, scale, displace_vector):
        self.scale, self.displace_vector = scale, displace_vector
        if self.si_point_data['coordinates'] is None:
            self.si_point_data['coordinates'] = (self.point_data['coordinates'] - self.displace_vector) / self.scale