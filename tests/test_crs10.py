import unittest
from VehicleData import build_falcon9
from vehicle_parameters import Falcon9Parameters


class CRS10Checks(unittest.TestCase):
    def parameters(self):
        return Falcon9Parameters()

    def test_budget_excludes_fairing_and_does_not_add_cargo_twice(self):
        parameters = self.parameters()
        vehicle = build_falcon9(parameters)
        self.assertFalse(hasattr(vehicle, "fairing_mass"))
        self.assertEqual(parameters.cargo_mass_kg, 2490)
        self.assertEqual(vehicle.getmass(), 553030)
        self.assertAlmostEqual(vehicle.get_mass_properties().mass_kg, vehicle.getmass())

    def test_separation_keeps_dragon_on_upper_vehicle(self):
        vehicle = build_falcon9(self.parameters())
        booster, upper = vehicle.separate()
        self.assertEqual(booster.payload_mass, 0)
        self.assertEqual(upper.payload_mass, 8430)
        self.assertEqual(booster.getmass(), 433100)
        self.assertEqual(booster.getmass() + upper.getmass(), vehicle.getmass())

    def test_default_cylinder_and_cg_use_the_same_geometry(self):
        parameters = self.parameters()
        self.assertEqual(parameters.mission, "CRS-10")
        self.assertAlmostEqual(parameters.payload_length, 7.2)
        self.assertAlmostEqual(parameters.payload_radius, 1.83)
        self.assertAlmostEqual(parameters.payload_cg_station, 58.8)
        self.assertAlmostEqual(parameters.height, 62.4)
        self.assertAlmostEqual(parameters.mass, build_falcon9().getmass())


if __name__ == "__main__":
    unittest.main()
