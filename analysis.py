import h5py
import numpy as np

# Path to your DREAM.3D file
file_path = 'data/forsterite.dream3d'

def print_structure(name, obj):
    # This function will be called for every item in the HDF5 file
    # `name` is the path to the item
    # `obj` is the actual object (group or dataset)
    if isinstance(obj, h5py.Dataset):
        print(f"Dataset: {name}, Shape: {obj.shape}, Type: {obj.dtype}")
    else:
        print(f"Group: {name}")

# Open the file
with h5py.File(file_path, 'r') as file:
    # Recursively list all groups and datasets in the file
    file.visititems(print_structure)

def extract_data(file_path):
    # Open the DREAM.3D file
    with h5py.File(file_path, 'r') as file:
        # Extract the EulerAngles, Phases, CrystalStructures, and LatticeConstants
        eulerAngles = file['DataContainers/ImageDataContainer/CellData/EulerAngles'][:]
        phases = file['DataContainers/ImageDataContainer/CellData/Phases'][:]
        crystalStructures = file['DataContainers/ImageDataContainer/CellEnsembleData/CrystalStructures'][:]
        latticeConstants = file['DataContainers/ImageDataContainer/CellEnsembleData/LatticeConstants'][:]

        return eulerAngles, phases, crystalStructures, latticeConstants


eulerAngles, phases, crystalStructures, latticeConstants = extract_data(file_path)

print("EulerAngles:", eulerAngles.shape)
print("Phases:", phases.shape)
print("CrystalStructures:", crystalStructures.shape)
print("LatticeConstants:", latticeConstants.shape)


# Saving data to a text file that can be read by Fortran
np.savetxt('abaqus/EulerAngles.txt', eulerAngles.reshape(-1, 3), header='phi1 PHI phi2', comments='')
np.savetxt('abaqus/Phases.txt', phases.flatten(), header='Phase index', comments='')
