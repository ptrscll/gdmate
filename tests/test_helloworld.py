"""
Test for placeholder code.
"""
import unittest

import gdmate.helloworld

class TestHelloWorld(unittest.TestCase):
    def test_helloworld(self):
        self.assertEqual(gdmate.helloworld.helloworld(),'Hello World')

if __name__ == '__main__':
    unittest.main()