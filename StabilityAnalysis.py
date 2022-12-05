import pymesh
import numpy as np
import matplotlib.pyplot as plt
import os
from Buoy import Buoy
import time
import csv
import imageio

class StabilityAnalysis():
    """Contains relevant methods to perform stability analysis on `Buoy` object"""
    
    def __init__(self, buoy:Buoy, angular_resolution:int, buoyancy_accuracy:float, rotation_axis:"list[float]", filename:str): # buoyancy accuracy, angular resolution, axis of rotation
        self.buoy = buoy
        self.angular_resolution = angular_resolution
        self.buoyancy_accuracy = buoyancy_accuracy
        self.rotation_axis = rotation_axis
        self.filename = filename

        self._zero()

    def _zero(self) -> None:
        """Translate `self.mesh` to be centered on X and Y axes, minimum Z value = 0."""
        mins, maxs = self.buoy.mesh.bbox[0], self.buoy.mesh.bbox[1]

        # first zero in global goordinate frame:
        offsets = [-(maxs[i] + mins[i])/2 for i in [0,1]] # center model on x and y axes
        offsets.append(-mins[2]) # move min z value to zero

        offsets = [offsets[i] - self.buoy.CG[i] for i in [0,1,2]] # move CG to origin
        self.buoy.mesh = self._translate(self.buoy.mesh, offsets) # translate to center
                
    def _translate(self, obj: pymesh.Mesh, distance:"list[float]") -> pymesh.Mesh:
        """Translate arbitrary mesh by `distance`. `distance` must be passed in form `[x,y,z]`. Return translated mesh."""
        return pymesh.form_mesh(obj.vertices + [distance], obj.faces, obj.voxels)

    def _rotate(self, obj: pymesh.Mesh, angle:int, vector) -> None:
        """Rotate mesh by `angle` about `vector`, rotating through the CG."""
        quat = pymesh.Quaternion.fromAxisAngle(vector, np.radians(angle))
        vertices = [quat.rotate(v) for v in obj.vertices]
        return pymesh.form_mesh(vertices, obj.faces, obj.voxels)

    def _findVolume(self, obj: pymesh.Mesh) -> float:
        """Return total volume of input `obj` mesh."""
        assert(obj.dim == 3)
        total_vol = 0.0

        for f in obj.faces:
            corners = obj.vertices[f]
            volume = np.dot(np.cross(corners[0], corners[1]), corners[2])
            total_vol += volume

        total_vol /= 6.0
        return total_vol

    def _findCentroid(self, obj: pymesh.Mesh) -> float:
        """Return centroid of input `obj` mesh."""
        assert(obj.dim == 3)
        volumes = np.zeros(obj.faces.shape[0])
        centroids = np.zeros([obj.faces.shape[0], 3])
        
        for i,f in enumerate(obj.faces):
            corners = obj.vertices[f]
            volumes[i] = np.dot(np.cross(corners[0], corners[1]), corners[2]) / 6
            centroids[i] = np.vstack([corners, np.zeros(3)]).mean(axis=0) # average all corners of tetrahedron, including the origin

        return (volumes * centroids.T).sum(axis=1) / (volumes.sum())

    def _sliceAtWaterline(self, obj: pymesh.Mesh) -> "tuple[pymesh.Mesh, float]":
        """Use `self.mass` to locate water line of `obj` mesh. Return sliced mesh and height of water line."""
        return self._sliceAtWaterlineRec(obj, obj.bbox[0][2], obj.bbox[1][2]) # use recursive function to find waterline

    def _sliceAtWaterlineRec(self, obj: pymesh.Mesh, lower_bound:float, upper_bound:float) -> "tuple[pymesh.Mesh, float]":
        """Recursively find the height of the waterline and return."""
        slice_height = (upper_bound+lower_bound)/2 # check the volume in the center of the two bounds
        mins = np.append(obj.bbox[0][0:2], slice_height) # min z value is slice_height
        maxs = obj.bbox[1][0:3]
        
        box = pymesh.generate_box_mesh(mins, maxs) # create box mesh encompassing area to remove
        sliced_obj = pymesh.boolean(obj, box, "difference", engine="igl") # subtract box from buoy mesh
        
        buoyancy = self._mm3_to_kg(self._findVolume(sliced_obj))
        if (np.abs(buoyancy - self.buoy.mass) < self.buoyancy_accuracy):
            return (sliced_obj, slice_height)
        elif (buoyancy < self.buoy.mass):
            return self._sliceAtWaterlineRec(obj, slice_height, upper_bound) # if volume is too small, search again in upper half
        else:
            return self._sliceAtWaterlineRec(obj, lower_bound, slice_height) # else volume is too big, so search again in lower half

    def _mm3_to_kg(self, volume:float):
        """Return mass of input volume (mm^3) of seawater"""
        return volume / 10**9 * 1020 # convert to m^3, then multiply by density of seawater

    def stabilityAnalysis(self) -> bool:
        """Perform full stability analysis of `mesh`, populating `self.CB` and `self.water_line`. Return True if buoy is expected to be stable."""
        if self._mm3_to_kg(self._findVolume(self.buoy.mesh)) < self.buoy.mass:
            print("Buoy will sink!")
            return False
        else:
            print("Beginning stability analysis")

        for angle in np.arange(0, 360, self.angular_resolution):
            t0 = time.time()
            rot_mesh = self._rotate(self.buoy.mesh, angle, self.rotation_axis)
            sliced_mesh, self.buoy.water_line[angle] = self._sliceAtWaterline(rot_mesh)
            self.buoy.CB[angle] = self._findCentroid(sliced_mesh)
            print(f"angle {angle} process time: {round(time.time() - t0, 2)} seconds")

        return self._checkStability()

    def _checkStability(self) -> bool:
        """Assess stability based on center of buoyancy position at each tilt angle"""
        axis = int(self.rotation_axis == [1,0,0]) # axis=1 (y axis) if rotating about x axis, otherwise axis=0 (x axis)
        moment_arms = {}
        self.moments = {}
        for angle in self.buoy.CB:
            if angle <= 180:
                moment_arms[angle] = self.buoy.CG[axis] - self.buoy.CB[angle][axis]
            else:
                moment_arms[angle] = self.buoy.CB[angle][axis] - self.buoy.CG[axis]

            self.moments[angle] = (moment_arms[angle] / 1000.0) * (self.buoy.mass * 9.81) # righting moment at each tilt angle

        ignore_instability = 5 # define angular region where unstable righting moments are tolerable
        for angle in self.buoy.CB:
            stable_region = (abs(angle%180) > ignore_instability) and (abs(angle%180 - 180) > ignore_instability)

            if stable_region and (moment_arms[angle] < 0):
                return False

        return True

    def plotRightingMoment(self):
        """Make a plot with angular tilt on the X-axis, and the righting moment on the Y-axis"""
        plt.plot(self.moments.keys(), self.moments.values(), marker='.', markersize=5)
        plt.xlabel("Angular Tilt (degrees)")
        plt.ylabel("Righting Moment (N-m)")
        plt.title("Buoyancy Righting Moment vs. Tilt Angle")
        plt.savefig(f"output/{self.filename}/righting_moments.jpg", dpi=300)
        plt.close()

    def writeToCSV(self):
        """Populates a CSV with the results of the stability analysis"""
        print("Writing results to CSV")
        fields = ["Angle (degrees)", "CB (mm)", "Water Line Height (mm)", "Righting Moment (N-m)"]
        rows = [[angle, self.buoy.CB[angle], self.buoy.water_line[angle], self.moments[angle]] for angle in self.buoy.CB]

        with open(f"output/{self.filename}/results.csv", 'w') as csvfile:
            csvwriter = csv.writer(csvfile) 
            csvwriter.writerow(fields) 
            csvwriter.writerows(rows)

    def visualize(self):
        """Make an animation showing buoy tilt and waterline"""
        max_dim = round(np.max(self.buoy.mesh.bbox)) * 1.5
        water = pymesh.generate_box_mesh([-max_dim]*3, [max_dim, max_dim, 0])

        print("Creating plots...")
        for angle in self.buoy.CB:
            # transform buoy mesh into correct position
            trans_mesh = self._translate(self._rotate(self.buoy.mesh, angle, self.rotation_axis), [0, 0, -self.buoy.water_line[angle]])

            # make plot:
            ax = plt.figure().gca(projection='3d')
            ax.set_box_aspect(aspect = (1,1,1))
            ax.set_position([0, 0, 1, 1])

            ax.plot_trisurf(water.vertices[:,0], water.vertices[:,1], water.vertices[:,2], triangles=water.faces, linewidth=0.2, antialiased=True, color='tab:blue', alpha=0.2)
            ax.plot_trisurf(trans_mesh.vertices[:,0], trans_mesh.vertices[:,1], trans_mesh.vertices[:,2], triangles=trans_mesh.faces, linewidth=0.2, color='grey', alpha=0.4, antialiased=True)
            ax.scatter(self.buoy.CG[0], self.buoy.CG[1], -self.buoy.water_line[angle], color='orangered', label="Center of Gravity")
            ax.scatter(self.buoy.CB[angle][0], self.buoy.CB[angle][1], self.buoy.CB[angle][2]-self.buoy.water_line[angle], color='royalblue', label="Center of Buoyancy")

            ax.set_xlim(-max_dim, max_dim)
            ax.set_ylim(-max_dim, max_dim)
            ax.set_zlim(-max_dim, max_dim)
            ax.set_xlabel("X (mm)")
            ax.set_ylabel("Y (mm)")
            ax.set_zlabel("Z (mm)")
            ax.legend(loc="upper right")
            
            # view plot along rotation axis
            if self.rotation_axis == [1,0,0]:
                ax.view_init(0, 0)
            else:
                ax.view_init(0, 270)

            plt.savefig(f"output/{self.filename}/frames/{angle}.jpg", dpi=400)
            plt.close()

        # Create GIF:
        print("Generating GIF...")
        image_files = [f"output/{self.filename}/frames/{angle}.jpg" for angle in self.buoy.CB]
        frames  = []
        frame_rate = len(self.buoy.CB) / 5.0

        for image_name in image_files:
            image = imageio.imread(image_name)
            frames.append(image)

        imageio.mimsave(f"output/{self.filename}/animation.gif", frames, 'GIF', fps=frame_rate)

        for image_name in image_files:
            os.remove(image_name)
