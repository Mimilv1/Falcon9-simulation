"""Team integration entry point. Flight events have not been selected yet."""
from VehicleData import build_falcon9


def main():
    vehicle = build_falcon9()
    print(f"Vehicle model loaded: {vehicle.getmass():,.0f} kg.")
    print("Flight simulation pending: connect team loads and agreed flight events.")


if __name__ == "__main__":
    main()
