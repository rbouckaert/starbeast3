package starbeast3.util;

import beast.base.util.Randomizer;

/**
 *  this class is used to make sure the apache library uses random numbers from the BEAST Randomizer
 *  so that the MCMC chain remains deterministic and starting with a certain seed twice will result in
 *  the same sequence.
 */
public class MyRandomizer implements org.apache.commons.math3.random.RandomGenerator {

	@Override
	public double nextDouble() {
		return Randomizer.nextDouble();
	}

	@Override
	public float nextFloat() {
		return Randomizer.nextFloat();
	}

	@Override
	public int nextInt() {
		return Randomizer.nextInt();
	}

	@Override
	public long nextLong() {
		return Randomizer.nextLong();
	}

	@Override
	public void setSeed(int seed) {
		Randomizer.setSeed(seed);			
	}

	@Override
	public void setSeed(int[] seed) {
		throw new RuntimeException("Not implemented");
		// Randomizer.setSeed(seed);			
	}

	@Override
	public void setSeed(long seed) {
		Randomizer.setSeed(seed);			
	}

	@Override
	public void nextBytes(byte[] bytes) {
		Randomizer.nextBytes(bytes);
	}

	@Override
	public int nextInt(int n) {
		return Randomizer.nextInt(n);
	}

	@Override
	public boolean nextBoolean() {
		return Randomizer.nextBoolean();
	}

	@Override
	public double nextGaussian() {
		return Randomizer.nextGaussian();
	}
	
}
