import unittest
from VehicleData import build_falcon9
from vehicle_parameters import CRS10Parameters


class CRS10Checks(unittest.TestCase):
    def parameters(self):
        return CRS10Parameters(dragon_loaded_mass_kg=9000, dragon_cg_station_m=58,
                               dragon_equivalent_length_m=6, dragon_equivalent_radius_m=1.83)

    def test_budget_excludes_fairing_and_does_not_add_cargo_twice(self):
        parameters = self.parameters()
        vehicle = build_falcon9(parameters)
        self.assertFalse(vehicle.fairing_attached)
        self.assertEqual(vehicle.fairing_mass, 0)
        self.assertEqual(parameters.cargo_mass_kg, 2490)
        self.assertEqual(vehicle.getmass(), 544600 + 9000)
        self.assertAlmostEqual(vehicle.get_mass_properties().mass_kg, vehicle.getmass())

    def test_separation_keeps_dragon_on_upper_vehicle(self):
        vehicle = build_falcon9(self.parameters())
        booster, upper = vehicle.separate()
        self.assertEqual(booster.payload_mass, 0)
        self.assertEqual(upper.payload_mass, 9000)
        self.assertEqual(booster.getmass(), 433100)
        self.assertEqual(booster.getmass() + upper.getmass(), vehicle.getmass())

    def test_unknown_mass_and_geometry_are_not_silently_assumed(self):
        with self.assertRaises(TypeError):
            CRS10Parameters()
        with self.assertRaises(ValueError):
            CRS10Parameters(dragon_loaded_mass_kg=2490, dragon_cg_station_m=58,
                             dragon_equivalent_length_m=6, dragon_equivalent_radius_m=1.83)


if __name__ == "__main__":
    unittest.main()
