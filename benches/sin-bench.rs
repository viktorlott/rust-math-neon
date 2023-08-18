use criterion::{black_box, criterion_group, criterion_main, Criterion};

const TEST_VALUES: &[f32] = &[
    -std::f32::consts::PI,
    -std::f32::consts::FRAC_PI_2,
    0.0,
    std::f32::consts::FRAC_PI_4,
    std::f32::consts::FRAC_PI_2,
    std::f32::consts::PI,
    2.0 * std::f32::consts::PI,
    0.4324545,
    0.26565545,
    0.060666,
    0.000434343,
    0.06066623,
    0.00043434343,
];
const VECTOR_SIZE: usize = 4; // Typically, NEON operations work on 4 single-precision floats at once.

fn bench_sinf_(c: &mut Criterion) {
    // Prepare a vector of test values.
    let test_values_vector: Vec<f32> = Vec::from(TEST_VALUES);

    // Ensure the size is a multiple of the vector size (e.g., 4 for NEON).
    let mut extended_test_values = test_values_vector.clone();
    while extended_test_values.len() % VECTOR_SIZE != 0 {
        extended_test_values.push(0.0); // Pad with zeroes.
    }

    let mut result_vector = vec![0.0_f32; extended_test_values.len()];

    c.bench_function("sinfv_neon", |b| {
        b.iter(|| {
            unsafe {
                rust_math_neon::sinfv_neon(
                    extended_test_values.as_mut_ptr(),
                    extended_test_values.len() as std::os::raw::c_int,
                    result_vector.as_mut_ptr(),
                );
            }
            black_box(&result_vector);
        })
    });
}

fn bench_rust_sin(c: &mut Criterion) {
    // Prepare a vector of test values.
    let test_values_vector: Vec<f32> = Vec::from(TEST_VALUES);

    let mut result_vector = vec![0.0_f32; test_values_vector.len()];

    c.bench_function("rust_sin", |b| {
        b.iter(|| {
            for (i, &val) in test_values_vector.iter().enumerate() {
                result_vector[i] = val.sin();
            }
            black_box(&result_vector);
        })
    });
}

criterion_group!(benches, bench_sinf_, bench_rust_sin);
criterion_main!(benches);
