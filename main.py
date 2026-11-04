"""Team integration entry point. Flight events have not been selected yet."""
import numpy as np
from frames import body_to_ecef_matrix
from VehicleData import build_falcon9


def main():
    vehicle = build_falcon9()
    state = vehicle.to_state()
    position_ecef = body_to_ecef_matrix(state.quaternion_body_to_ecef_wxyz) @ state.position_body_m
    print(f"{vehicle.parameters.mission} model loaded: {vehicle.getmass():,.0f} kg.")
    print("ECEF convention: Earth centre at (0,0,0); launch site on the +z axis.")
    print("Initial state: base at (0, 0, {:,.0f}) m; CG at ({:,.1f}, {:,.1f}, {:,.1f}) m.".format(
        vehicle.position[0], *(np.round(position_ecef, 6) + 0.0)))
    print("Flight simulation pending: connect team loads and agreed flight events.")


if __name__ == "__main__":
    main()
