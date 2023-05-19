declare module 'ml-airpls' {
  export interface AirPLSOptions {
    /**
     * @default 100
     */
    maxIterations?: number;
    /**
     * @default 100
     */
    lambda?: number;
    /**
     * @default 0.001
     */
    factorCriterion?: number;
    /**
     * array of number to weight the baseline. by default all point has the same weight.
     */
    weights?: number[],
    /**
     * array of indexes of points that should have the maximum weight on each iteration.
     * @default []
     */
    controlPoints?: number[],
    /**
     * array of indexes of points that should have the maximum weight on each iteration.
     * @default []
     */
    baseLineZones?: Array<{from: number, to: number}>,
  }

  export default function airPLS(
    x: number[] | Float64Array,
    y: number[] | Float64Array,
    options?: AirPLSOptions,
  ): {
    corrected: number[] | Float64Array,
    baseline: number[] | Float64Array,
    iteration: number,
    error: number
  };
}