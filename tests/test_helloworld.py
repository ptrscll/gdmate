"""
Test for placeholder code.
"""
import unittest

import gdmate as gd

class TestHelloWorld(unittest.TestCase):
    def test_helloworld(self):
        self.assertEqual(gd.helloworld(),'Hello World')

if __name__ == '__main__':
    unittest.main()