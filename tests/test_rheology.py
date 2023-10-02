"""
Test for rheology.py
"""
import unittest

import gdmate as gd

class TestHelloWorld(unittest.TestCase):
    def test_cond_geotherm(self):
        self.assertEqual(gd.helloworld(),'Hello World')
    
    def test_adiab_geotherm(self):
        self.assertEqual(gd.helloworld(),'Hello World')
    
    def test_geotherm(self):
        self.assertEqual(gd.helloworld(),'Hello World')


if __name__ == '__main__':
    unittest.main()