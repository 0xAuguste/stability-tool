from Buoy import Buoy
from StabilityAnalysis import StabilityAnalysis
import argparse
import os

if __name__ == "__main__":
    # Process command line arguments:
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--filename', type=str, required=True, help="Name of the mesh file")
    parser.add_argument('-m', '--mass', type=float, required=True, help="Mass of the buoy in kilograms.")
    parser.add_argument('-c', '--cg', type=float, nargs="+", required=True, help="Center of gravity location (x y z).")
    parser.add_argument('-r', '--resolution', type=int, default=45, help="Angular resolution (in degrees) of stability analysis.")
    parser.add_argument('-a', '--accuracy', type=float, default=0.001, help="Buoyancy accuracy of water line computation (in kilograms). Default is 1 gram accuracy, i.e. buoyancy will match mass within 1 gram.")
    parser.add_argument('-v', '--rotation', type=str, default='x', help="Rotational axis to spin the buoy about. Options are 'x' or 'y'.")
    args = parser.parse_args()

    # Interpret rotational axis:
    if args.rotation.lower() == 'x':
        axis = [1,0,0]
    elif args.rotation.lower() == 'y':
        axis = [0,-1,0]
    else:
        raise ValueError("Invalid rotational axis.")

    # Create necessary output directory:
    filename = os.path.splitext(args.filename)[0]
    if not os.path.exists(f'output/{filename}'):
        os.makedirs(f'output/{filename}')
    os.makedirs(f'output/{filename}/frames')

    # Run the analysis:
    buoy = Buoy(f"input/{args.filename}", args.mass, args.cg)
    sim = StabilityAnalysis(buoy, args.resolution, args.accuracy, axis, filename)
    try:
        stable = sim.stabilityAnalysis()
    except:
        os.removedirs(f'output/{filename}/frames')
    if stable:
        print("Buoy is stable")
    else:
        print("Buoy is unstable")
    sim.plotRightingMoment()
    sim.writeToCSV()
    sim.visualize()
    os.removedirs(f'output/{filename}/frames')
    print("Done.")