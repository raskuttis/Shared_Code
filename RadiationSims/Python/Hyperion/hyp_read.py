import subprocess
import re
import sys
import numpy as np
import string as string
import struct
from ath_units import *
from scipy.interpolate import interp1d

class AthenaDataRemote():
    def __init__(self,hostname,datafolder,basename,tnum,Nproc,units=None):
        
        fname = '{0:s}id{2:d}/{1:s}.{3:04d}.vtk'.format(datafolder,basename,0,tnum)
        self.filenames = [fname]
        self.Nproc = Nproc
        for step in range(1,Nproc):
            fname = '{0:s}id{2:d}/{1:s}-id{2:d}.{3:04d}.vtk'.format(datafolder,basename,step,tnum)
            self.filenames.append(fname)
        self.hostname = hostname

        fname = self.filenames[0]
        sshout = subprocess.Popen(['ssh', hostname, 'cat', fname],
                                  stdout=subprocess.PIPE)
                                  
        alllines = sshout.stdout.readlines()
        nlines = np.size(alllines)
        nline = 0
        line = alllines[nline]
        grid = {}
        if units is not None:
            self.units = units
        else:
            # Initialize physical units:
            self.units = Units_LMT()

        while line != '':
            splitup = line.strip().split()
            _parse_vtk(splitup, grid)
            if "CELL_DATA" in splitup:
                break
            nline = nline + 1
            line = alllines[nline]


        self.coordsys = 'Cartesian'
        self.gamma = 5.0/3.0
        self.iso_csound = 2.0
        self.mu = 2.33*self.units.mH
            
        self.time  = np.array(grid['time']).astype(float)
        self.Nx    = np.array(grid['dimensions']).astype(int) - 1
        self.ndim  = np.array(sum(self.Nx > 0)).astype(int)
        self.Npx = [int(round(Nproc**(1.0/self.ndim))) for x in range(0,self.ndim)]
        self.Nx = self.Nx * self.Npx
        self.dx    = np.array(grid['dx'])
        self.MinX  = np.array(grid['origin'])
        self.MaxX  = self.MinX + self.Nx*self.dx
        self.x1nodes = np.linspace(self.MinX[0],self.MaxX[0],self.Nx[0]+1)
        self.x2nodes = np.linspace(self.MinX[1],self.MaxX[1],self.Nx[1]+1)
        self.x3nodes = np.linspace(self.MinX[2],self.MaxX[2],self.Nx[2]+1)
        self.x1zones = np.linspace(self.MinX[0]+0.5*self.dx[0],self.MaxX[0]-0.5*self.dx[0],self.Nx[0])
        self.x2zones = np.linspace(self.MinX[1]+0.5*self.dx[1],self.MaxX[1]-0.5*self.dx[1],self.Nx[1])
        self.x3zones = np.linspace(self.MinX[2]+0.5*self.dx[2],self.MaxX[2]-0.5*self.dx[2],self.Nx[2])
            
        self.ncells = np.array(grid['ncells']).astype(int) * self.Nproc
        if self.Nx[self.Nx > 0].prod() != self.ncells:
            raise ValueError('Number of cells = {}, but Nx = {}'.format(self.ncells,self.Nx))

        self._field_map = {}
        nlinestart = nline
        nline = nline + 1
        line = alllines[nline]
        while line != '' and nline < nlines - 1:
            splitup = line.strip().split()
            _parse_vtk(splitup, grid)
            if "SCALARS" in splitup:
                field = grid['read_field']
                nline = nline + 1
                line = alllines[nline]
                self._field_map[field] = ('scalar', nline - 1 - nlinestart)
            elif "VECTORS" in splitup:
                field = splitup[1]
                self._field_map[field] = ('vector', nline - nlinestart)
            nline = nline + 1
            line = alllines[nline]

        self.fields = self._field_map.keys()
        self.mhd = 'cell_centered_B' in self.fields
        self.adiabatic = 'total_energy' in self.fields
        self.self_gravity = 'gravitational_potential' in self.fields
        self.radiation = 'rad_energy_density' in self.fields

    def get_field(self,field='density'):
        
        field_type, linenum = self._field_map[field]
        if   field_type=='scalar':
            nvar = 1
        elif field_type=='vector':
            nvar = 3
        datadims = np.append(self.Nx,nvar)
        data = np.zeros(datadims)
        if field in ('density','momentum','total_energy','rad_energy_density',
                     'rad_flux','gravitational_potential'):
            Nproc = self.Nproc
            for step in range(0,Nproc):
                print step
                hostname = self.hostname
                filename = self.filenames[step]
                sshout = subprocess.Popen(['ssh', hostname, 'cat', filename],
                                          stdout=subprocess.PIPE)
                grid = {}
                nline = 0
                line = sshout.stdout.readline()
                while line != '':
                    splitup = line.strip().split()
                    _parse_vtk(splitup, grid)
                    if "CELL_DATA" in splitup:
                        break
                    nline = nline + 1
                    line = sshout.stdout.readline()
                locNx = np.array(grid['dimensions']).astype(int) - 1
                locncells = np.array(grid['ncells']).astype(int)
                locMinX = np.array(grid['origin'])
                locivals = self._xyz_to_ijk(xcut = locMinX[0], ycut = locMinX[1], zcut = locMinX[2])
                locMinijk = [locivals['icut'],locivals['jcut'],locivals['kcut']]
                locMaxijk = locMinijk + locNx
                
                for i in range(0,linenum):
                    line = sshout.stdout.readline()

                splitup = line.strip().split()
                _parse_vtk(splitup, grid)
                #print splitup
                if "SCALARS" in splitup:
                    fieldread = grid['read_field']
                    nvar = 1
                    line = sshout.stdout.readline()  # Read the lookup table line
                elif "VECTORS" in splitup:
                    fieldread = splitup[1]
                    nvar = 3
                else:
                    fieldread = ' '
                    nvar = 0
                #print fieldread, field
                if fieldread==field:
                    fmt = '>{}f'.format(nvar*locncells)
                    size = struct.calcsize(fmt)
                else:
                    print 'Error in reading data'
                    exit()

                locdata = sshout.stdout.read(size)
                dims = locNx
                dims[locNx==0] = 1
                dims = np.append(np.flipud(dims),nvar)
                locdata = np.transpose(np.reshape(np.array(struct.unpack(fmt,locdata)),dims),(2,1,0,3))
                data[locMinijk[0]:locMaxijk[0],locMinijk[1]:locMaxijk[1],locMinijk[2]:locMaxijk[2],:] = locdata
                
            if nvar==1:
                data = np.squeeze(data,3)
            return data
                     
        elif field=='momentum_x':
            return self.get_field('momentum')[:,:,:,0]
        elif field=='momentum_y':
            return self.get_field('momentum')[:,:,:,1]
        elif field=='momentum_z':
            return self.get_field('momentum')[:,:,:,2]
        elif field=='cell_centered_Bx':
            return self.get_field('cell_centered_B')[:,:,:,0]
        elif field=='cell_centered_By':
            return self.get_field('cell_centered_B')[:,:,:,1]
        elif field=='cell_centered_Bz':
            return self.get_field('cell_centered_B')[:,:,:,2]
        elif field=='rad_flux_x':
            return self.get_field('rad_flux')[:,:,:,0]
        elif field=='rad_flux_y':
            return self.get_field('rad_flux')[:,:,:,1]
        elif field=='rad_flux_z':
            return self.get_field('rad_flux')[:,:,:,2]
        elif field=='velocity_x':
            return self.get_field('momentum_x')/self.get_field('density')
        elif field=='velocity_y':
            return self.get_field('momentum_y')/self.get_field('density')
        elif field=='velocity_z':
            return self.get_field('momentum_z')/self.get_field('density')
        elif field=='velocity':
            return self.get_field('momentum')/np.expand_dims(self.get_field('rho'),axis=3)
        elif field=='momentum_magnitude':
            return _magnitude(self.get_field('momentum'))
        elif field=='velocity_magnitude':
            return _magnitude(self.get_field('velocity'))
        elif field=='cell_centered_B_magnitude':
            return _magnitude(self.get_field('cell_centered_B'))
        elif field=='rad_flux_magnitude':
            return _magnitude(self.get_field('rad_flux'))
        elif field=='kinetic_energy':
            rho = self.get_field('density')
            return 0.5*np.sum(self.get_field('momentum')**2,axis=3)/rho
        elif field=='internal_energy':
            E = self.get_field('total_energy')
            Emag = self.get_field('magnetic_pressure')
            Ekin = self.get_field('kinetic_energy')
            return E - Emag - Ekin
        elif field=='pressure':
            return self.get_field('internal_energy')*(self.gamma - 1.0)
        elif field=='magnetic_pressure':
            if self.mhd:
                return 0.5*np.sum(self.get_field('cell_centered_B')**2,axis=3)
            else:
                return 0.0
        elif field=='temperature':
            P = self.get_field('pressure')
            rho = self.get_field('density')
            return (self.mu/self.units.kB)*P/rho
        elif field=='Mach_number':
            return self.get_field('velocity_magnitude')/self.get_field('sound_speed')
        elif field=='sound_speed':
            if self.adiabatic:
                P = self.get_field('pressure')
                rho = self.get_field('density')
                return np.sqrt(self.gamma*P/rho)
            else:
                return self.iso_csound
        elif field=='reduced_flux_magnitude':
            F = self.get_field('rad_flux_magnitude')
            Er = self.get_field('rad_energy_density')
            return F/(self.units.c*Er)
        elif field=='rad_temperature':
            return (self.get_field('rad_energy_density')/self.units.aR)**0.25
        elif field=='Er':
            return self.get_field('rad_energy_density')
        elif field=='e':
            return self.get_field('internal_energy')
        elif field in ('F','Frad'):
            return self.get_field('rad_flux_magnitude')
        elif field=='B':
            return self.get_field('cell_centered_B_magnitude')
        elif field in ('P','Pgas','Pg'):
            return self.get_field('pressure')
        elif field in ('v','V'):
            return self.get_field('velocity_magnitude')
        elif field in ('vx','Vx','v1','V1'):
            return self.get_field('velocity_x')
        elif field in ('vy','Vy','v2','V2'):
            return self.get_field('velocity_y')
        elif field in ('vz','Vz','v3','V3'):
            return self.get_field('velocity_z')
        elif field in ('Fx','F1'):
            return self.get_field('rad_flux_x')
        elif field in ('Fy','F2'):
            return self.get_field('rad_flux_y')
        elif field in ('Fz','F3'):
            return self.get_field('rad_flux_z')
        elif field in ('Bx','B1'):
            return self.get_field('cell_centered_Bx')
        elif field in ('By','B2'):
            return self.get_field('cell_centered_By')
        elif field in ('Bz','B3'):
            return self.get_field('cell_centered_Bz')
        elif field in ('rho','d'):
            return self.get_field('density')
        elif field=='cs':
            return self.get_field('sound_speed')
        elif field=='f':
            return self.get_field('reduced_flux_magnitude')
        elif field in ('Tgas','T','Tg'):
            return self.get_field('temperature')
        elif field in ('Trad','Tr'):
            return self.get_field('rad_temperature')
        elif field=='Tdiff':
            return self.get_field('temperature')-self.get_field('rad_temperature')
        else:
            raise ValueError('Unknown field {}'.format(field))
                
    def get_slice_ijk(self,field='density',transpose=False,**kwargs):
        if   'icut' in kwargs.keys():
            imask = kwargs['icut']
            jmask = slice(0,self.Nx[1])
            kmask = slice(0,self.Nx[2])
            x,y = self.x2nodes,self.x3nodes
        elif 'jcut' in kwargs.keys():
            imask = slice(0,self.Nx[0])
            jmask = kwargs['jcut']
            kmask = slice(0,self.Nx[2])
            x,y = self.x1nodes,self.x3nodes
        elif 'kcut' in kwargs.keys():
            imask = slice(0,self.Nx[0])
            jmask = slice(0,self.Nx[1])
            kmask = kwargs['kcut']
            x,y = self.x1nodes,self.x2nodes
        else:
            raise ValueError('Exactly one of icut,jcut,kcut must be defined!')
        Z = self.get_field(field)[imask,jmask,kmask]
        if transpose:
            X,Y = np.meshgrid(y,x,indexing='ij')
            Z = Z.T
        else:
            X,Y = np.meshgrid(x,y,indexing='ij')
        return X,Y,Z

        
    def get_lineout_ijk(self,field='density',**kwargs):
        if   all(x in kwargs.keys() for x in ['jcut','kcut']):
            imask = slice(0,self.Nx[0])
            jmask = kwargs['jcut']
            kmask = kwargs['kcut']
            x = self.x1zones
        elif all(x in kwargs.keys() for x in ['icut','kcut']):
            imask = kwargs['icut']
            jmask = slice(0,self.Nx[1])
            kmask = kwargs['kcut']
            x = self.x2zones
        elif all(x in kwargs.keys() for x in ['icut','jcut']):
            imask = kwargs['icut']
            jmask = kwargs['jcut']
            kmask = slice(0,self.Nx[2])
            x = self.x3zones
        else:
            raise ValueError('Exactly two of icut,jcut,kcut must be defined!')
        y = self.get_field(field)[imask,jmask,kmask]
        return x,y

    def get_slice_xyz(self,field='density',transpose=False,**kwargs):
        kwargs = self._xyz_to_ijk(**kwargs)
        return self.get_slice_ijk(field,transpose,**kwargs)
        
    def get_lineout_xyz(self,field='density',**kwargs):
        kwargs = self._xyz_to_ijk(**kwargs)
        return self.get_lineout_ijk(field,**kwargs)

    def _xyz_to_ijk(self,**kwargs):
        kwargs_out = {}
        if 'xcut' in kwargs.keys():
            if self.dx[0]>0.0:
                icut = np.floor_divide(kwargs['xcut']-self.MinX[0],self.dx[0]).astype(int)
            else:
                icut = 0
            kwargs_out['icut'] = icut
        if 'ycut' in kwargs.keys():
            if self.dx[1]>0.0:
                jcut = np.floor_divide(kwargs['ycut']-self.MinX[1],self.dx[1]).astype(int)
            else:
                jcut = 0
            kwargs_out['jcut'] = jcut
        if 'zcut' in kwargs.keys():
            if self.dx[2]>0.0:
                kcut = np.floor_divide(kwargs['zcut']-self.MinX[2],self.dx[2]).astype(int)
            else:
                kcut = 0        
            kwargs_out['kcut'] = kcut
        return kwargs_out

class AthenaDataRemoteJoined():
    def __init__(self,hostname,datafolder,basename,tnum,units=None):
        
        fname = '{0:s}{1:s}_joined.{2:04d}.vtk'.format(datafolder,basename,tnum)
        self.filename = fname
        self.hostname = hostname
        sshout = subprocess.Popen(['ssh', hostname, 'cat', fname],
                                  stdout=subprocess.PIPE)
                                  
        grid = {}
        if units is not None:
            self.units = units
        else:
            # Initialize physical units:
            self.units = Units_LMT()

        line = sshout.stdout.readline()
        while line != '':
            splitup = line.strip().split()
            _parse_vtk(splitup, grid)
            if "CELL_DATA" in splitup:
                break
            line = sshout.stdout.readline()

        self.coordsys = 'Cartesian'
        self.gamma = 5.0/3.0
        self.iso_csound = 2.0
        self.mu = 2.33*self.units.mH
            
        self.time  = np.array(grid['time']).astype(float)
        self.Nx    = np.array(grid['dimensions']).astype(int) - 1
        self.ndim  = np.array(sum(self.Nx > 0)).astype(int)
        self.dx    = np.array(grid['dx'])
        self.MinX  = np.array(grid['origin'])
        self.MaxX  = self.MinX + self.Nx*self.dx
        self.x1nodes = np.linspace(self.MinX[0],self.MaxX[0],self.Nx[0]+1)
        self.x2nodes = np.linspace(self.MinX[1],self.MaxX[1],self.Nx[1]+1)
        self.x3nodes = np.linspace(self.MinX[2],self.MaxX[2],self.Nx[2]+1)
        self.x1zones = np.linspace(self.MinX[0]+0.5*self.dx[0],self.MaxX[0]-0.5*self.dx[0],self.Nx[0])
        self.x2zones = np.linspace(self.MinX[1]+0.5*self.dx[1],self.MaxX[1]-0.5*self.dx[1],self.Nx[1])
        self.x3zones = np.linspace(self.MinX[2]+0.5*self.dx[2],self.MaxX[2]-0.5*self.dx[2],self.Nx[2])
        
        self.ncells = np.array(grid['ncells']).astype(int)
        if self.Nx[self.Nx > 0].prod() != self.ncells:
            raise ValueError('Number of cells = {}, but Nx = {}'.format(self.ncells,self.Nx))

    def get_field(self,field='density'):
        
        if field in ('density','momentum','total_energy','rad_energy_density',
                     'rad_flux','gravitational_potential'):
            hostname = self.hostname
            filename = self.filename
            sshout = subprocess.Popen(['ssh', hostname, 'cat', filename],
                                        stdout=subprocess.PIPE)
            
            grid = {}
            line = sshout.stdout.readline()
            nline = 0
            while line != '':
                #print nline
                nline = nline + 1
                splitup = line.strip().split()
                _parse_vtk(splitup, grid)
                if "SCALARS" in splitup:
                    fieldread = grid['read_field']
                    nvar = 1
                    line = sshout.stdout.readline()  # Read the lookup table line
                elif "VECTORS" in splitup:
                    fieldread = splitup[1]
                    nvar = 3
                else:
                    fieldread = ' '
                    nvar = 0
                if fieldread==field:
                    fmt = '>{}f'.format(nvar*self.ncells)
                    size = struct.calcsize(fmt)
                    break
                else:
                    line = sshout.stdout.readline()

            data = sshout.stdout.read(size)
            dims = self.Nx
            dims[self.Nx==0] = 1
            dims = np.append(np.flipud(dims),nvar)
            data = np.transpose(np.reshape(np.array(struct.unpack(fmt,data)),dims),(2,1,0,3))
                
            if nvar==1:
                data = np.squeeze(data,3)
            return data
                     
        elif field=='momentum_x':
            return self.get_field('momentum')[:,:,:,0]
        elif field=='momentum_y':
            return self.get_field('momentum')[:,:,:,1]
        elif field=='momentum_z':
            return self.get_field('momentum')[:,:,:,2]
        elif field=='cell_centered_Bx':
            return self.get_field('cell_centered_B')[:,:,:,0]
        elif field=='cell_centered_By':
            return self.get_field('cell_centered_B')[:,:,:,1]
        elif field=='cell_centered_Bz':
            return self.get_field('cell_centered_B')[:,:,:,2]
        elif field=='rad_flux_x':
            return self.get_field('rad_flux')[:,:,:,0]
        elif field=='rad_flux_y':
            return self.get_field('rad_flux')[:,:,:,1]
        elif field=='rad_flux_z':
            return self.get_field('rad_flux')[:,:,:,2]
        elif field=='velocity_x':
            return self.get_field('momentum_x')/self.get_field('density')
        elif field=='velocity_y':
            return self.get_field('momentum_y')/self.get_field('density')
        elif field=='velocity_z':
            return self.get_field('momentum_z')/self.get_field('density')
        elif field=='velocity':
            return self.get_field('momentum')/np.expand_dims(self.get_field('rho'),axis=3)
        elif field=='momentum_magnitude':
            return _magnitude(self.get_field('momentum'))
        elif field=='velocity_magnitude':
            return _magnitude(self.get_field('velocity'))
        elif field=='cell_centered_B_magnitude':
            return _magnitude(self.get_field('cell_centered_B'))
        elif field=='rad_flux_magnitude':
            return _magnitude(self.get_field('rad_flux'))
        elif field=='kinetic_energy':
            rho = self.get_field('density')
            return 0.5*np.sum(self.get_field('momentum')**2,axis=3)/rho
        elif field=='internal_energy':
            E = self.get_field('total_energy')
            Emag = self.get_field('magnetic_pressure')
            Ekin = self.get_field('kinetic_energy')
            return E - Emag - Ekin
        elif field=='pressure':
            return self.get_field('internal_energy')*(self.gamma - 1.0)
        elif field=='magnetic_pressure':
            if self.mhd:
                return 0.5*np.sum(self.get_field('cell_centered_B')**2,axis=3)
            else:
                return 0.0
        elif field=='temperature':
            P = self.get_field('pressure')
            rho = self.get_field('density')
            return (self.mu/self.units.kB)*P/rho
        elif field=='Mach_number':
            return self.get_field('velocity_magnitude')/self.get_field('sound_speed')
        elif field=='sound_speed':
            if self.adiabatic:
                P = self.get_field('pressure')
                rho = self.get_field('density')
                return np.sqrt(self.gamma*P/rho)
            else:
                return self.iso_csound
        elif field=='reduced_flux_magnitude':
            F = self.get_field('rad_flux_magnitude')
            Er = self.get_field('rad_energy_density')
            return F/(self.units.c*Er)
        elif field=='rad_temperature':
            return (self.get_field('rad_energy_density')/self.units.aR)**0.25
        elif field=='Er':
            return self.get_field('rad_energy_density')
        elif field=='e':
            return self.get_field('internal_energy')
        elif field in ('F','Frad'):
            return self.get_field('rad_flux_magnitude')
        elif field=='B':
            return self.get_field('cell_centered_B_magnitude')
        elif field in ('P','Pgas','Pg'):
            return self.get_field('pressure')
        elif field in ('v','V'):
            return self.get_field('velocity_magnitude')
        elif field in ('vx','Vx','v1','V1'):
            return self.get_field('velocity_x')
        elif field in ('vy','Vy','v2','V2'):
            return self.get_field('velocity_y')
        elif field in ('vz','Vz','v3','V3'):
            return self.get_field('velocity_z')
        elif field in ('Fx','F1'):
            return self.get_field('rad_flux_x')
        elif field in ('Fy','F2'):
            return self.get_field('rad_flux_y')
        elif field in ('Fz','F3'):
            return self.get_field('rad_flux_z')
        elif field in ('Bx','B1'):
            return self.get_field('cell_centered_Bx')
        elif field in ('By','B2'):
            return self.get_field('cell_centered_By')
        elif field in ('Bz','B3'):
            return self.get_field('cell_centered_Bz')
        elif field in ('rho','d'):
            return self.get_field('density')
        elif field=='cs':
            return self.get_field('sound_speed')
        elif field=='f':
            return self.get_field('reduced_flux_magnitude')
        elif field in ('Tgas','T','Tg'):
            return self.get_field('temperature')
        elif field in ('Trad','Tr'):
            return self.get_field('rad_temperature')
        elif field=='Tdiff':
            return self.get_field('temperature')-self.get_field('rad_temperature')
        else:
            raise ValueError('Unknown field {}'.format(field))
                
    def get_slice_ijk(self,field='density',transpose=False,**kwargs):
        if   'icut' in kwargs.keys():
            imask = kwargs['icut']
            jmask = slice(0,self.Nx[1])
            kmask = slice(0,self.Nx[2])
            x,y = self.x2nodes,self.x3nodes
        elif 'jcut' in kwargs.keys():
            imask = slice(0,self.Nx[0])
            jmask = kwargs['jcut']
            kmask = slice(0,self.Nx[2])
            x,y = self.x1nodes,self.x3nodes
        elif 'kcut' in kwargs.keys():
            imask = slice(0,self.Nx[0])
            jmask = slice(0,self.Nx[1])
            kmask = kwargs['kcut']
            x,y = self.x1nodes,self.x2nodes
        else:
            raise ValueError('Exactly one of icut,jcut,kcut must be defined!')
        Z = self.get_field(field)[imask,jmask,kmask]
        if transpose:
            X,Y = np.meshgrid(y,x,indexing='ij')
            Z = Z.T
        else:
            X,Y = np.meshgrid(x,y,indexing='ij')
        return X,Y,Z

        
    def get_lineout_ijk(self,field='density',**kwargs):
        if   all(x in kwargs.keys() for x in ['jcut','kcut']):
            imask = slice(0,self.Nx[0])
            jmask = kwargs['jcut']
            kmask = kwargs['kcut']
            x = self.x1zones
        elif all(x in kwargs.keys() for x in ['icut','kcut']):
            imask = kwargs['icut']
            jmask = slice(0,self.Nx[1])
            kmask = kwargs['kcut']
            x = self.x2zones
        elif all(x in kwargs.keys() for x in ['icut','jcut']):
            imask = kwargs['icut']
            jmask = kwargs['jcut']
            kmask = slice(0,self.Nx[2])
            x = self.x3zones
        else:
            raise ValueError('Exactly two of icut,jcut,kcut must be defined!')
        y = self.get_field(field)[imask,jmask,kmask]
        return x,y

    def get_slice_xyz(self,field='density',transpose=False,**kwargs):
        kwargs = self._xyz_to_ijk(**kwargs)
        return self.get_slice_ijk(field,transpose,**kwargs)
        
    def get_lineout_xyz(self,field='density',**kwargs):
        kwargs = self._xyz_to_ijk(**kwargs)
        return self.get_lineout_ijk(field,**kwargs)

    def _xyz_to_ijk(self,**kwargs):
        kwargs_out = {}
        if 'xcut' in kwargs.keys():
            if self.dx[0]>0.0:
                icut = np.floor_divide(kwargs['xcut']-self.MinX[0],self.dx[0]).astype(int)
            else:
                icut = 0
            kwargs_out['icut'] = icut
        if 'ycut' in kwargs.keys():
            if self.dx[1]>0.0:
                jcut = np.floor_divide(kwargs['ycut']-self.MinX[1],self.dx[1]).astype(int)
            else:
                jcut = 0
            kwargs_out['jcut'] = jcut
        if 'zcut' in kwargs.keys():
            if self.dx[2]>0.0:
                kcut = np.floor_divide(kwargs['zcut']-self.MinX[2],self.dx[2]).astype(int)
            else:
                kcut = 0        
            kwargs_out['kcut'] = kcut
        return kwargs_out


class AthenaData():
    def __init__(self,vtkfile):
        self.units = Units_LMT()
        self.filename = vtkfile
        
        with open(self.filename,'rb') as f:
            grid = {}
            line = f.readline()
            while line != '':
                splitup = line.strip().split()
                _parse_vtk(splitup, grid)
                if "CELL_DATA" in splitup:
                    break
                line = f.readline()
            
            # These 3 parameters can be overwritten...
            self.coordsys = 'Cartesian'
            self.gamma = 5.0/3.0
            self.iso_csound = 2.0
            # So can these
            self.mu = 2.33*self.units.mH
            
            self.time  = np.array(grid['time']).astype(float)
            self.Nx    = np.array(grid['dimensions']).astype(int) - 1
            self.dx    = np.array(grid['dx'])
            self.MinX  = np.array(grid['origin'])
            self.MaxX  = self.MinX + self.Nx*self.dx
            self.ndim  = np.array(sum(self.Nx > 0)).astype(int)
            self.x1nodes = np.linspace(self.MinX[0],self.MaxX[0],self.Nx[0]+1)
            self.x2nodes = np.linspace(self.MinX[1],self.MaxX[1],self.Nx[1]+1)
            self.x3nodes = np.linspace(self.MinX[2],self.MaxX[2],self.Nx[2]+1)
            self.x1zones = np.linspace(self.MinX[0]+0.5*self.dx[0],self.MaxX[0]-0.5*self.dx[0],self.Nx[0])
            self.x2zones = np.linspace(self.MinX[1]+0.5*self.dx[1],self.MaxX[1]-0.5*self.dx[1],self.Nx[1])
            self.x3zones = np.linspace(self.MinX[2]+0.5*self.dx[2],self.MaxX[2]-0.5*self.dx[2],self.Nx[2])
            
            self.ncells = np.array(grid['ncells']).astype(int)
            if self.Nx[self.Nx > 0].prod() != self.ncells:
                raise ValueError('Number of cells = {}, but Nx = {}'.format(self.ncells,self.Nx))
            
            self._field_map = {}
            data_offset = f.tell()  # This is just after the "CELL_DATA" line
            line = f.readline()
            while line != '':
                splitup = line.strip().split()
                _parse_vtk(splitup, grid)
                if "SCALARS" in splitup:
                    field = grid['read_field']
                    line = f.readline()  # Read the lookup table line
                    self._field_map[field] = ('scalar',f.tell()-data_offset)
                elif "VECTORS" in splitup:
                    field = splitup[1]
                    self._field_map[field] = ('vector',f.tell()-data_offset)
                line = f.readline()

        self.fields = self._field_map.keys()
        self.mhd = 'cell_centered_B' in self.fields
        self.adiabatic = 'total_energy' in self.fields
        self.self_gravity = 'gravitational_potential' in self.fields
        self.radiation = 'rad_energy_density' in self.fields

    def get_field(self,field='density'):
        if field in ('density','momentum','total_energy','rad_energy_density',
                     'rad_flux','gravitational_potential'):
            with open(self.filename,'rb') as f:
                line = f.readline()
                while line != '':
                    splitup = line.strip().split()
                    if "CELL_DATA" in splitup:
                        break
                    line = f.readline()
                data_offset = f.tell()
                field_type,field_offset = self._field_map[field]
                f.seek(data_offset+field_offset)
                if   field_type=='scalar':
                    nvar = 1
                elif field_type=='vector':
                    nvar = 3
                fmt = '>{}f'.format(nvar*self.ncells)
                size = struct.calcsize(fmt)
                data = f.read(size)
            dims = self.Nx
            dims[self.Nx==0] = 1
            dims = np.append(np.flipud(dims),nvar)
            data = np.transpose(np.reshape(np.array(struct.unpack(fmt,data)),dims),(2,1,0,3))
            if nvar==1:
                data = np.squeeze(data,3)
            return data
        elif field=='momentum_x':
            return self.get_field('momentum')[:,:,:,0]
        elif field=='momentum_y':
            return self.get_field('momentum')[:,:,:,1]
        elif field=='momentum_z':
            return self.get_field('momentum')[:,:,:,2]
        elif field=='cell_centered_Bx':
            return self.get_field('cell_centered_B')[:,:,:,0]
        elif field=='cell_centered_By':
            return self.get_field('cell_centered_B')[:,:,:,1]
        elif field=='cell_centered_Bz':
            return self.get_field('cell_centered_B')[:,:,:,2]
        elif field=='rad_flux_x':
            return self.get_field('rad_flux')[:,:,:,0]
        elif field=='rad_flux_y':
            return self.get_field('rad_flux')[:,:,:,1]
        elif field=='rad_flux_z':
            return self.get_field('rad_flux')[:,:,:,2]
        elif field=='velocity_x':
            return self.get_field('momentum_x')/self.get_field('density')
        elif field=='velocity_y':
            return self.get_field('momentum_y')/self.get_field('density')
        elif field=='velocity_z':
            return self.get_field('momentum_z')/self.get_field('density')
        elif field=='velocity':
            return self.get_field('momentum')/np.expand_dims(self.get_field('rho'),axis=3)
        elif field=='momentum_magnitude':
            return _magnitude(self.get_field('momentum'))
        elif field=='velocity_magnitude':
            return _magnitude(self.get_field('velocity'))
        elif field=='cell_centered_B_magnitude':
            return _magnitude(self.get_field('cell_centered_B'))
        elif field=='rad_flux_magnitude':
            return _magnitude(self.get_field('rad_flux'))
        elif field=='kinetic_energy':
            rho = self.get_field('density')
            return 0.5*np.sum(self.get_field('momentum')**2,axis=3)/rho
        elif field=='internal_energy':
            E = self.get_field('total_energy')
            Emag = self.get_field('magnetic_pressure')
            Ekin = self.get_field('kinetic_energy')
            return E - Emag - Ekin
        elif field=='pressure':
            return self.get_field('internal_energy')*(self.gamma - 1.0)
        elif field=='magnetic_pressure':
            if self.mhd:
                return 0.5*np.sum(self.get_field('cell_centered_B')**2,axis=3)
            else:
                return 0.0
        elif field=='temperature':
            P = self.get_field('pressure')
            rho = self.get_field('density')
            return (self.mu/self.units.kB)*P/rho
        elif field=='Mach_number':
            return self.get_field('velocity_magnitude')/self.get_field('sound_speed')
        elif field=='sound_speed':
            if self.adiabatic:
                P = self.get_field('pressure')
                rho = self.get_field('density')
                return np.sqrt(self.gamma*P/rho)
            else:
                return self.iso_csound
        elif field=='reduced_flux_magnitude':
            F = self.get_field('rad_flux_magnitude')
            Er = self.get_field('rad_energy_density')
            return F/(self.units.c*Er)
        elif field=='rad_temperature':
            return (self.get_field('rad_energy_density')/self.units.aR)**0.25
        elif field=='Er':
            return self.get_field('rad_energy_density')
        elif field=='e':
            return self.get_field('internal_energy')
        elif field in ('F','Frad'):
            return self.get_field('rad_flux_magnitude')
        elif field=='B':
            return self.get_field('cell_centered_B_magnitude')
        elif field in ('P','Pgas','Pg'):
            return self.get_field('pressure')
        elif field in ('v','V'):
            return self.get_field('velocity_magnitude')
        elif field in ('vx','Vx','v1','V1'):
            return self.get_field('velocity_x')
        elif field in ('vy','Vy','v2','V2'):
            return self.get_field('velocity_y')
        elif field in ('vz','Vz','v3','V3'):
            return self.get_field('velocity_z')
        elif field in ('Fx','F1'):
            return self.get_field('rad_flux_x')
        elif field in ('Fy','F2'):
            return self.get_field('rad_flux_y')
        elif field in ('Fz','F3'):
            return self.get_field('rad_flux_z')
        elif field in ('Bx','B1'):
            return self.get_field('cell_centered_Bx')
        elif field in ('By','B2'):
            return self.get_field('cell_centered_By')
        elif field in ('Bz','B3'):
            return self.get_field('cell_centered_Bz')
        elif field in ('rho','d'):
            return self.get_field('density')
        elif field=='cs':
            return self.get_field('sound_speed')
        elif field=='f':
            return self.get_field('reduced_flux_magnitude')
        elif field in ('Tgas','T','Tg'):
            return self.get_field('temperature')
        elif field in ('Trad','Tr'):
            return self.get_field('rad_temperature')
        elif field=='Tdiff':
            return self.get_field('temperature')-self.get_field('rad_temperature')
        else:
            raise ValueError('Unknown field {}'.format(field))

def _parse_vtk(splitup,grid):
    # grid is a dictionary
    if "vtk" in splitup:
        grid['vtk_version'] = splitup[-1]
    elif "time=" in splitup:
        time_index = splitup.index("time=")
        grid['time'] = float(splitup[time_index+1].rstrip(','))
        grid['level'] = int(splitup[time_index+3].rstrip(','))
        grid['domain'] = int(splitup[time_index+5].rstrip(','))
    elif "DIMENSIONS" in splitup:
        grid['dimensions'] = np.array(splitup[-3:]).astype('int')
    elif "ORIGIN" in splitup:
        grid['origin'] = np.array(splitup[-3:]).astype('float64')
    elif "SPACING" in splitup:
        grid['dx'] = np.array(splitup[-3:]).astype('float64')
    elif "CELL_DATA" in splitup:
        grid["ncells"] = int(splitup[-1])
    elif "SCALARS" in splitup:
        field = splitup[1]
        grid['read_field'] = field
        grid['read_type'] = 'scalar'
    elif "VECTORS" in splitup:
        field = splitup[1]
        grid['read_field'] = field
        grid['read_type'] = 'vector'

def read_fluxfile(hostname,fluxfile):
    sshflux = subprocess.Popen(['ssh', hostname, 'cat', fluxfile],
                              stdout=subprocess.PIPE)
    data = np.loadtxt(sshflux.stdout)
    dsize = np.shape(data)
    nvs = dsize[0]
    nfs = dsize[1]
    tvals = data[:,0]
    uniqt = np.unique(tvals)
    nts = len(uniqt)
    nrs = nvs / nts
    outdata = np.reshape(data, (nts, nrs, nfs))
    return outdata

def read_outfile(hostname,outfile):
    sshout = subprocess.Popen(['ssh', hostname, 'cat', outfile],
                              stdout=subprocess.PIPE)
    alllines = sshout.stdout.readlines()
    return alllines

def read_local_outfile(outfile):
    f = open(outfile, 'r')
    alllines = f.readlines()
    f.close()
    return alllines

def read_hstfile(hostname,hstfile):
    sshhst = subprocess.Popen(['ssh', hostname, 'cat', hstfile],
                              stdout=subprocess.PIPE)
    data = np.loadtxt(sshhst.stdout)
    return data

def read_denfloorfile(hostname,denflfile):
    sshfloor = subprocess.Popen(['ssh', hostname, 'cat', denflfile],
                              stdout=subprocess.PIPE)
    data = np.loadtxt(sshfloor.stdout)
    return data

def read_local_hstfile(hstfile):
    f = open(hstfile, 'r')
    data = np.loadtxt(f)
    f.close()
    return data

def read_pdffile(hostname,outfile):
    sshout = subprocess.Popen(['ssh', hostname, 'cat', outfile],
                              stdout=subprocess.PIPE)
    alllines = sshout.stdout.readlines()
    return alllines

def search_stardir(hostname,starbase):
    lfile = 1
    nstar = 0
    while (lfile > 0):
        starfile = starbase + string.zfill(str(nstar),4) + '.dat'
        lsout = subprocess.Popen(['ssh', hostname, 'cat', starfile],
                                 stdout=subprocess.PIPE)
        lfile = len(lsout.stdout.readlines())
        nstar = nstar + 1

    return nstar - 1

def read_allstars(hostname,starbase,time,tff):
    nstars = search_stardir(hostname,starbase)
    nts = len(time)
    allstars = np.zeros((nts,nstars,11))
    for i in xrange(1,nstars+1):
        starfile = starbase + string.zfill(str(i-1),4) + '.dat'
        data = read_starfile(hostname,starfile)
        ndims = len(np.shape(data))
        if ndims > 1:
            tdat = data[:,0] / tff
        else:
            tdat = data / tff
        allstars[:,i-1,0] = time
        for j in xrange(1,11):
            if ndims > 1:
                xdat = data[:,j]
                xf = interp1d(tdat, xdat, bounds_error=False, fill_value=0)
                xtime = xf(time)
                allstars[:,i-1,j] = xtime
            else:
                xdata = data

    return allstars

def read_starfile(hostname,starfile):
    lsout = subprocess.Popen(['ssh', hostname, 'cat', starfile],
                             stdout=subprocess.PIPE)
    data = np.loadtxt(lsout.stdout)
    return data

def search_local_stardir(starbase):
    lfile = 1
    nstar = 0
    while (lfile > 0):
        starfile = starbase + string.zfill(str(nstar),4) + '.dat'
        try:
            lsout = open(starfile,'r')
            lfile = len(lsout.readlines())
            lsout.close()
        except:
            lfile = 0
        nstar = nstar + 1

    return nstar - 1

def read_local_allstars(starbase,time,tff):
    nstars = search_local_stardir(starbase)
    nts = len(time)
    allstars = np.zeros((nts,nstars,11))
    for i in xrange(1,nstars+1):
        starfile = starbase + string.zfill(str(i-1),4) + '.dat'
        data = read_local_starfile(starfile)
        ndims = len(np.shape(data))
        if ndims > 1:
            tdat = data[:,0] / tff
        else:
            tdat = data / tff
        allstars[:,i-1,0] = time
        for j in xrange(1,11):
            if ndims > 1:
                xdat = data[:,j]
                xf = interp1d(tdat, xdat, bounds_error=False, fill_value=0)
                xtime = xf(time)
                allstars[:,i-1,j] = xtime
            else:
                xdata = data

    return allstars

def read_local_starfile(starfile):
    lsout = open(starfile,'r')
    data = np.loadtxt(lsout)
    lsout.close()
    return data

def float2explatex(val):
    exp = np.log10(val)
    exp = int(exp)
    mant = int(np.round(val / 10**exp))
    strout = '${0:d} \\times 10^{1:d}$'.format(mant,exp)
    return strout

def modelstr(mcloud, rcloud):
    exp = np.log10(mcloud)
    exp = int(exp+0.05)
    mant = int(np.round(mcloud / 10**exp))
    rout = int(rcloud)
    strout = '$\\Sigma$-M{0:d}E{1:d}-R{2:d}'.format(mant,exp,rout)
    return strout