import pymesh
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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

    def _rotate(self, obj: pymesh.Mesh, angle:int, vector) -> pymesh.Mesh:
        """Rotate mesh by `angle` about `vector`, rotating through the origin."""
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
            volumes[i] = np.dot(np.cross(corners[0], corners[1]), corners[2]) / 6.0
            centroids[i] = np.vstack([corners, np.zeros(3)]).mean(axis=0) # average all corners of tetrahedron, including the origin

        return (volumes * centroids.T).sum(axis=1) / (volumes.sum()) # weighted average of all centroids

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

    def _mm3_to_kg(self, volume:float) -> float:
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
            moment_arms[angle] = self.buoy.CB[angle][axis]
            self.moments[angle] = (moment_arms[angle] / 1000.0) * (self.buoy.mass * 9.81) # righting moment at each tilt angle

        ignore_instability = 5 # define angular region where unstable righting moments are tolerable
        for angle in self.buoy.CB:
            stable_region = (abs(angle%180) > ignore_instability) and (abs(angle%180 - 180) > ignore_instability)

            if stable_region and ((0 <= angle <= 180 and moment_arms[angle] > 0) or (180 < angle <= 360 and moment_arms[angle] < 0)) :
                return False

        return True

    def findRideAngle(self) -> "list[float]":
        """Returns the angle the buoy will sit in static water. If multiple stable angles exist, the angle closest to vertical is returned."""
        stable_zero_crossings = []
        stable_angles = []
        angles = list(self.moments.keys())
        moments = list(self.moments.values())
        for i in range(len(angles)):
            if (round(moments[i], 4) < 0) and (round(moments[i-1], 4) > 0):
                stable_zero_crossings.append([angles[i-1], angles[i]])
            elif (round(moments[i], 4) == 0) and (round(moments[i-1], 4) > 0) and (round(moments[i+1], 4) < 0):
                stable_angles.append(angles[i])

        for crossing_angles in stable_zero_crossings:
            angles = crossing_angles.copy()
            if crossing_angles[0] > crossing_angles[1]:
                angles[0] -=360

            # Interpolate between points to find precise zero crossing:
            diff = -self.moments[crossing_angles[0]] * (angles[1] - angles[0]) / (self.moments[crossing_angles[1]] - self.moments[crossing_angles[0]])
            stable_angles.append(angles[0] + diff)

        ride_angles = [angle if angle <= 180 else angle-360 for angle in stable_angles] # map range [0,360] to [-180,180]
        return ride_angles

    def plotRightingMoment(self) -> None:
        """Make plots to visualize the righting moment"""
        angles = np.array(list(self.moments.keys()))
        moments = np.array([round(moment, 4) for moment in self.moments.values()])
        angles = np.append(angles, 360)
        moments = np.append(moments, moments[0])
        colors =  [1 if moment < 0 else -0.25 if moment == 0 else -1 for moment in moments] # create color map based on righting moment

        self._plotCartesianRightingMoment(angles, moments, colors)
        self._plotPolarRightingMoment(angles, moments, colors)

    def _plotCartesianRightingMoment(self, angles:float, moments:float, colors:float):
        """Make a plot with angular tilt on the X-axis, and the righting moment on the Y-axis"""
        cmap = plt.cm.Spectral
        legend_points = [Line2D([0], [0], color=cmap(1.), lw=4), Line2D([0], [0], color=cmap(-1.), lw=4)]
        plt.plot(angles, moments, color="darkgrey", linewidth=0.5, zorder=1)
        plt.scatter(angles, moments, marker='.', c=colors, cmap="Spectral", zorder=2)
        plt.xlabel("Angular Tilt (degrees)")
        plt.ylabel("Moment (N-m)")
        plt.title("Buoyancy Moment About CG vs. Tilt Angle")
        plt.legend(legend_points, ["CW moment", "CCW moment"])
        plt.savefig(f"output/{self.filename}/righting_moments.jpg", dpi=300)
        plt.close()

    def _plotPolarRightingMoment(self, angles:float, moments:float, colors:float):
        """Make a plot with angular tilt on the theta axis, and the righting moment on the radius axis"""
        cmap = plt.cm.Spectral
        legend_points = [Line2D([0], [0], color=cmap(1.), lw=4), Line2D([0], [0], color=cmap(-1.), lw=4)]
        radians = angles*np.pi/180
        fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
        ax.plot(radians, moments, color="darkgrey", linewidth=0.5, zorder=1)
        ax.scatter(radians, moments, marker='.', c=colors, cmap="Spectral", zorder=2)
        ax.grid(True)
        ax.set_theta_zero_location('S')
        ax.set_title("Buoyancy Moment (N-m) About CG")
        ax.legend(legend_points, ["CW moment", "CCW moment"], bbox_to_anchor=(-.3, -0.1), loc="lower left")
        plt.savefig(f"output/{self.filename}/righting_moments_polar.jpg", dpi=300)
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
            ax.scatter(0, 0, -self.buoy.water_line[angle], color='orangered', label="Center of Gravity")
            ax.scatter(self.buoy.CB[angle][0], self.buoy.CB[angle][1], self.buoy.CB[angle][2]-self.buoy.water_line[angle], color='royalblue', label="Center of Buoyancy")

            ax.set_xlim(-max_dim, max_dim)
            ax.set_ylim(-max_dim, max_dim)
            ax.set_zlim(-max_dim, max_dim)
            ax.set_zlabel("Z (mm)")
            ax.legend(loc="upper right")
            
            # view plot along rotation axis
            if self.rotation_axis == [1,0,0]: # if rotating about x axis
                ax.view_init(0, 0)
                ax.set_xticks([])
                ax.set_ylabel("Y (mm)")
            else: # rotating about y axis
                ax.view_init(0, 270)
                ax.set_yticks([])
                ax.set_xlabel("X (mm)")

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
