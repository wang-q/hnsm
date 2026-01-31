// Add these imports to use the stdsimd library
#![feature(portable_simd)]
use std::simd::prelude::*;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use hnsm::{hash_hv_bit, hash_hv_i8};
use rand::rngs::SmallRng;
use rand::rngs::StdRng;
use rand::{Rng, RngCore, SeedableRng};
use rapidhash::{RapidHashSet, RapidRng};

pub fn encode_hash_hd_simd(seed_vec: &[u64], hv_d: usize) -> Vec<i16> {
    let num_seed = seed_vec.len();
    let mut hv = vec![-(num_seed as i16); hv_d];

    let num_chunk = hv_d / 64;

    // Loop through all seeds
    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(*hash);

        // SIMD-based HV encoding
        for i in 0..num_chunk {
            let rnd_bits = rng.next_u64();

            // Use SIMD to process 4 bits at a time
            for j in (0..64).step_by(4) {
                // Create a SIMD vector of 4 bits
                let bit_mask = u64x4::splat(1);
                let shift =
                    Simd::from_array([j as u64, (j + 1) as u64, (j + 2) as u64, (j + 3) as u64]);
                let bits = (u64x4::splat(rnd_bits) >> shift) & bit_mask;

                // Convert bits to i16 and shift left by 1
                let bits_i16 = bits.cast::<i16>() << Simd::splat(1);

                // Load the target HV values
                let mut hv_simd = i16x4::from_slice(&hv[i * 64 + j..i * 64 + j + 4]);

                // Accumulate the bits
                hv_simd += bits_i16;

                // Store the updated HV values
                hv_simd.copy_to_slice(&mut hv[i * 64 + j..i * 64 + j + 4]);
            }
        }
    }

    hv
}

pub fn encode_hash_hd_simd3(seed_vec: &[u64], hv_d: usize) -> Vec<i32> {
    let num_seed = seed_vec.len();
    let mut hv = vec![-(num_seed as i32); hv_d];

    let num_chunk = hv_d / 64;

    // Loop through all seeds
    for hash in seed_vec.iter() {
        let mut rng = RapidRng::seed_from_u64(*hash);

        // SIMD-based HV encoding
        for i in 0..num_chunk {
            let rnd_bits = rng.next_u64();
            let bits_ary: [u8; 8] = rnd_bits.to_be_bytes();

            // Use SIMD to process 8 bits at a time
            for (bi, j) in (0..64).step_by(8).enumerate() {
                let bit_mask = u8x8::splat(1);
                let shift = Simd::from_array([
                    j as u8,
                    (j + 1) as u8,
                    (j + 2) as u8,
                    (j + 3) as u8,
                    (j + 4) as u8,
                    (j + 5) as u8,
                    (j + 6) as u8,
                    (j + 7) as u8,
                ]);
                let bits = (u8x8::splat(bits_ary[bi]) >> shift) & bit_mask;
                let bits_i32 = bits.cast::<i32>() << Simd::splat(1);

                // Load the target HV values
                let mut hv_simd = i32x8::from_slice(&hv[i * 64 + j..i * 64 + j + 8]);

                // Accumulate the bits
                hv_simd += bits_i32;

                // Store the updated HV values
                hv_simd.copy_to_slice(&mut hv[i * 64 + j..i * 64 + j + 8]);
            }
        }
    }

    hv
}

pub fn encode_hash_hd_rapid(seed_vec: &[u64], hv_d: usize) -> Vec<i16> {
    let mut hv = vec![-(seed_vec.len() as i16); hv_d];

    for hash in seed_vec {
        let mut rng = RapidRng::seed_from_u64(*hash);

        for i in 0..(hv_d / 64) {
            let rnd_bits = rng.next_u64();

            for j in 0..64 {
                hv[i * 64 + j] += (((rnd_bits >> j) & 1) << 1) as i16;
            }
        }
    }

    hv
}

pub fn encode_hash_hd_std(seed_vec: &[u64], hv_d: usize) -> Vec<i16> {
    let mut hv = vec![-(seed_vec.len() as i16); hv_d];

    for hash in seed_vec {
        let mut rng = StdRng::seed_from_u64(*hash);

        for i in 0..(hv_d / 64) {
            let rnd_bits = rng.next_u64();

            for j in 0..64 {
                hv[i * 64 + j] += (((rnd_bits >> j) & 1) << 1) as i16;
            }
        }
    }

    hv
}

pub fn encode_hash_hd_small(seed_vec: &[u64], hv_d: usize) -> Vec<i16> {
    let mut hv = vec![-(seed_vec.len() as i16); hv_d];

    for hash in seed_vec {
        let mut rng = SmallRng::seed_from_u64(*hash);

        for i in 0..(hv_d / 64) {
            let rnd_bits = rng.next_u64();

            for j in 0..64 {
                hv[i * 64 + j] += (((rnd_bits >> j) & 1) << 1) as i16;
            }
        }
    }

    hv
}

// Generate a random k-mer hash set
fn generate_kmer_hash_set(size: usize) -> RapidHashSet<u64> {
    let mut rng = StdRng::seed_from_u64(42); // Fixed seed for reproducibility
    let mut kmer_hash_set = RapidHashSet::default();

    for _ in 0..size {
        kmer_hash_set.insert(rng.random::<u64>());
    }

    kmer_hash_set
}

// Benchmark function
fn bench_encode_hash_hd(c: &mut Criterion) {
    // Create test datasets of different sizes
    let kmer_hash_set_small = generate_kmer_hash_set(1000); // Small dataset
    let kmer_hash_set_medium = generate_kmer_hash_set(10_000); // Medium dataset

    let seed_vec_small: Vec<u64> = kmer_hash_set_small.iter().cloned().collect();
    let seed_vec_medium: Vec<u64> = kmer_hash_set_medium.iter().cloned().collect();

    let hv_d = 4096; // Hypervector dimension

    // Benchmark small dataset
    c.bench_function("encode_hash_hd_lib_small", |b| {
        b.iter(|| hash_hv_bit(black_box(&seed_vec_small), hv_d))
    });
    c.bench_function("encode_hash_hd_simd_small", |b| {
        b.iter(|| encode_hash_hd_simd(black_box(&seed_vec_small), hv_d))
    });
    c.bench_function("encode_hash_hd_simd_i8_small", |b| {
        b.iter(|| hash_hv_i8(black_box(&seed_vec_small), hv_d))
    });
    c.bench_function("encode_hash_hd_simd3_small", |b| {
        b.iter(|| encode_hash_hd_simd3(black_box(&seed_vec_small), hv_d))
    });
    c.bench_function("encode_hash_hd_rapid_small", |b| {
        b.iter(|| encode_hash_hd_rapid(black_box(&seed_vec_small), hv_d))
    });
    c.bench_function("encode_hash_hd_std_small", |b| {
        b.iter(|| encode_hash_hd_std(black_box(&seed_vec_small), hv_d))
    });
    c.bench_function("encode_hash_hd_small_small", |b| {
        b.iter(|| encode_hash_hd_small(black_box(&seed_vec_small), hv_d))
    });

    // Benchmark medium dataset
    c.bench_function("encode_hash_hd_lib_medium", |b| {
        b.iter(|| hash_hv_bit(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_simd_medium", |b| {
        b.iter(|| encode_hash_hd_simd(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_simd_i8_medium", |b| {
        b.iter(|| hash_hv_i8(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_simd3_medium", |b| {
        b.iter(|| encode_hash_hd_simd3(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_rapid_medium", |b| {
        b.iter(|| encode_hash_hd_rapid(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_std_medium", |b| {
        b.iter(|| encode_hash_hd_std(black_box(&seed_vec_medium), hv_d))
    });
    c.bench_function("encode_hash_hd_small_medium", |b| {
        b.iter(|| encode_hash_hd_small(black_box(&seed_vec_medium), hv_d))
    });
}

// Define benchmark group
criterion_group!(
    name = benches;
    config = Criterion::default().sample_size(10); // Set sample size
    targets = bench_encode_hash_hd
);

// Run benchmarks
criterion_main!(benches);

// encode_hash_hd_simd_small
//                         time:   [1.6922 ms 1.6984 ms 1.7052 ms]
//                         change: [-1.1949% -0.7620% -0.3383%] (p = 0.01 < 0.05)
//                         Change within noise threshold.

// encode_hash_hd_simd2_small
//                         time:   [676.41 µs 678.52 µs 679.93 µs]
//                         change: [+0.3887% +0.7028% +1.0266%] (p = 0.00 < 0.05)
//                         Change within noise threshold.

// encode_hash_hd_simd_i8_small
//                         time:   [419.74 µs 421.03 µs 422.67 µs]
// Found 1 outliers among 10 measurements (10.00%)
//   1 (10.00%) high mild

// encode_hash_hd_simd3_small
//                         time:   [1.5759 ms 1.5842 ms 1.5911 ms]
//                         change: [-1.9303% -0.7676% +0.3719%] (p = 0.25 > 0.05)
//                         No change in performance detected.
// Found 1 outliers among 10 measurements (10.00%)
//   1 (10.00%) high mild

// encode_hash_hd_rapid_small
//                         time:   [1.6542 ms 1.6630 ms 1.6775 ms]
//                         change: [+1.8280% +2.5541% +3.2866%] (p = 0.00 < 0.05)
//                         Performance has regressed.

// encode_hash_hd_std_small
//                         time:   [1.7369 ms 1.7491 ms 1.7643 ms]
//                         change: [+0.8124% +1.4188% +2.0306%] (p = 0.00 < 0.05)
//                         Change within noise threshold.

// encode_hash_hd_small_small
//                         time:   [1.6276 ms 1.6389 ms 1.6455 ms]
//                         change: [-2.8960% -1.6557% -0.4737%] (p = 0.02 < 0.05)
//                         Change within noise threshold.

// encode_hash_hd_simd_medium
//                         time:   [16.837 ms 16.937 ms 17.059 ms]
//                         change: [-0.4475% +0.1193% +0.7442%] (p = 0.71 > 0.05)
//                         No change in performance detected.

// encode_hash_hd_simd2_medium
//                         time:   [6.7508 ms 6.8296 ms 6.9440 ms]
//                         change: [-0.3764% +0.7010% +2.0167%] (p = 0.33 > 0.05)
//                         No change in performance detected.

// encode_hash_hd_simd_i8_medium
//                         time:   [4.1850 ms 4.1992 ms 4.2100 ms]

// encode_hash_hd_simd3_medium
//                         time:   [15.680 ms 15.820 ms 15.941 ms]
//                         change: [-3.0446% -1.3741% -0.1414%] (p = 0.09 > 0.05)
//                         No change in performance detected.

// encode_hash_hd_rapid_medium
//                         time:   [16.288 ms 16.415 ms 16.519 ms]
//                         change: [-1.3015% -0.4403% +0.4122%] (p = 0.34 > 0.05)
//                         No change in performance detected.
// Found 2 outliers among 10 measurements (20.00%)
//   2 (20.00%) high mild

// encode_hash_hd_std_medium
//                         time:   [17.393 ms 17.457 ms 17.514 ms]
//                         change: [-0.1103% +0.6447% +1.4475%] (p = 0.14 > 0.05)
//                         No change in performance detected.
// Found 1 outliers among 10 measurements (10.00%)
//   1 (10.00%) high mild

// encode_hash_hd_small_medium
//                         time:   [16.278 ms 16.338 ms 16.402 ms]
//                         change: [-1.4538% -0.6406% +0.1704%] (p = 0.16 > 0.05)
//                         No change in performance detected.

/*
# 基准测试分析总结 (2026-01-30) - 更新版

## 1. 变更说明
- **优化输入类型**: 所有函数（包括库函数 `hash_hv`）现在直接接收 `&[u64]` 向量，移除了函数内部从 `HashSet` 到 `Vec` 的转换开销。这回应了关于“冗余操作”的反馈。
- **基准测试调整**: `HashSet` 到 `Vec` 的转换现在在基准测试循环外进行，确保只测试核心编码算法的性能。

## 2. 性能排名 (从快到慢)

### 小数据集 (Small Dataset, 1,000 items)
1. **encode_hash_hd_simd_i8_small**: ~416 µs (新的 i8 实现)
2. encode_hash_hd_simd2_small: ~670 µs (当前库实现逻辑)
3. encode_hash_hd_lib_small: ~670 µs (直接调用库函数)
4. encode_hash_hd_simd3_small: ~1.53 ms
5. encode_hash_hd_small_small: ~1.65 ms
6. encode_hash_hd_rapid_small: ~1.66 ms
7. encode_hash_hd_simd_small: ~1.69 ms
8. encode_hash_hd_std_small: ~1.71 ms

### 中等数据集 (Medium Dataset, 10,000 items)
1. **encode_hash_hd_simd_i8_medium**: ~4.20 ms (新的 i8 实现)
2. encode_hash_hd_lib_medium: ~6.73 ms (直接调用库函数)
3. encode_hash_hd_simd2_medium: ~6.73 ms (当前库实现逻辑)
4. encode_hash_hd_simd3_medium: ~15.40 ms
5. encode_hash_hd_rapid_medium: ~16.50 ms
6. encode_hash_hd_small_medium: ~16.54 ms
7. encode_hash_hd_simd_medium: ~16.88 ms
8. encode_hash_hd_std_medium: ~17.16 ms

## 3. 主要观察

- **i8 方案保持领先**: 即使移除了输入转换的开销，`encode_hash_hd_simd_i8` 仍然是最快的，比当前的库实现 (`simd2`/`lib`) 快约 **1.6倍**。
- **库函数性能**: `encode_hash_hd_lib` 与本地优化的 `encode_hash_hd_simd2` 性能一致，证明库函数的开销极小。
- **消除冗余的影响**: 移除 `HashSet` -> `Vec` 的转换并没有改变各算法的相对排名，说明主要瓶颈在于 RNG 生成和位操作/SIMD 累加，而不是集合迭代。

## 4. 建议

- 推荐将 `encode_hash_hd_simd_i8` 的逻辑（使用 i8 累加代替位操作）合并到 `src/libs/hv.rs` 的主库中，以获得最大的性能提升。
- 保持 `hash_hv` 接收 `&[u64]` 的接口更改，为调用者提供更大的灵活性。
*/
