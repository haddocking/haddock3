"""Mathematical functions."""
import random


class RandomNumberGenerator:
    """
    Generate random numbers.

    Python uses the Mersenne Twister as the core generator.
    It produces 53-bit precision floats and has a period of 2**19937-1
    """

    def __init__(self, seed: int = 494) -> None:
        """494 = sum([ord(i) for i in 'HADDOCK'])."""
        self.seed = seed
        self.random = random.Random()
        self.random.seed(self.seed)

    def __call__(self,
                 lower_limit: float = 0.,
                 upper_limit: float = 1.) -> float:
        """Generate a random number between range."""
        return self.random.uniform(lower_limit, upper_limit)

    def randint(self, lower_limit: int = 0, upper_limit: int = 9) -> int:
        """Generate a random integer."""
        return int(self() * (upper_limit + 1)) + lower_limit
