use core::arch::aarch64::*;
use std::f32::consts::PI;

pub unsafe fn fft(inp: &[f32]) -> Vec<f32> {
    let n = inp.len();
    let zero = 0.0f32;

    if n == 1 {
        return vec![inp[0], zero];
    }
    // Instead of padding, we just naively run the dft.
    if n % 2 == 1 {
        return dft(inp);
    }

    let mut out = vec![zero; n * 2];
    let mut even = Vec::with_capacity(n / 2);
    let mut odd = Vec::with_capacity(n / 2);

    // Split the input into even and odd components.
    for (i, &inp) in inp.iter().enumerate() {
        if i % 2 == 0 {
            even.push(inp)
        } else {
            odd.push(inp);
        }
    }

    // Recursive FFT calls for even and odd parts.
    let even_fft = fft(&even);
    let odd_fft = fft(&odd);

    // Combine the FFT results from the even and odd parts.
    let two_pi = PI + PI;
    let n_t = n as f32;
    let n_half = n / 2;

    if n_half % 4 == 0 {
        vectorized_combine(two_pi, n_t, n_half, &odd_fft, &even_fft, &mut out);
    } else {
        for k in 0..n_half {
            let k_t = k as f32;
            let theta = two_pi * k_t / n_t;
            let re = theta.cos();
            let im = -theta.sin();

            let re_odd = odd_fft[2 * k];
            let im_odd = odd_fft[2 * k + 1];

            out[2 * k] = even_fft[2 * k] + re * re_odd - im * im_odd;
            out[2 * k + 1] = even_fft[2 * k + 1] + re * im_odd + im * re_odd;

            out[2 * (k + n_half)] = even_fft[2 * k] - re * re_odd + im * im_odd;
            out[2 * (k + n_half) + 1] = even_fft[2 * k + 1] - re * im_odd - im * re_odd;
        }
    }
    out
}

unsafe fn vectorized_combine(
    two_pi: f32,
    n_t: f32,
    n_half: usize,
    odd_fft: &Vec<f32>,
    even_fft: &Vec<f32>,
    out: &mut Vec<f32>,
) {
    let pi_vec = vdupq_n_f32(PI);
    let two_pi_vec = vdupq_n_f32(two_pi);
    let n_f32_vec = vdupq_n_f32(n_t);
    let increment = vld1q_f32([0.0, 1.0, 2.0, 3.0].as_ptr());

    for k in (0..n_half).step_by(4) {
        let k_base = vdupq_n_f32(k as f32);
        let k_values = vaddq_f32(k_base, increment);
        // pi * [k0..k3]
        let k_pi_values = vmulq_f32(two_pi_vec, k_values);

        // Compute theta_values using vectorized operations
        let s_theta = vdivq_f32(k_pi_values, n_f32_vec);

        // This is a trick such that we can calc cos with sin function. sin(x+pi/2) = cos(x)
        let cos_theta = vaddq_f32(s_theta, vdivq_f32(pi_vec, vdupq_n_f32(2.0)));

        let mut theta_output = [0.0f32; 8];
        vst1q_f32(theta_output.as_mut_ptr(), s_theta);
        vst1q_f32(theta_output.as_mut_ptr().add(4), cos_theta);

        let mut vector = [0.0f32; 8];
        rust_math_neon::sinfv_neon(theta_output.as_mut_ptr(), 8, vector.as_mut_ptr());

        let re_v = vld1q_f32(vector[4..].as_ptr());
        let im_v = vnegq_f32(vld1q_f32(vector[..4].as_ptr()));

        let re_x_odd = [
            odd_fft[2 * k],
            odd_fft[2 * (k + 1)],
            odd_fft[2 * (k + 2)],
            odd_fft[2 * (k + 3)],
        ];
        let im_x_odd = [
            odd_fft[2 * k + 1],
            odd_fft[2 * (k + 1) + 1],
            odd_fft[2 * (k + 2) + 1],
            odd_fft[2 * (k + 3) + 1],
        ];

        let re_x_even = [
            even_fft[2 * k],
            even_fft[2 * (k + 1)],
            even_fft[2 * (k + 2)],
            even_fft[2 * (k + 3)],
        ];
        let im_x_even = [
            even_fft[2 * k + 1],
            even_fft[2 * (k + 1) + 1],
            even_fft[2 * (k + 2) + 1],
            even_fft[2 * (k + 3) + 1],
        ];

        let re_odd_v = vld1q_dup_f32(re_x_odd.as_ptr());
        let im_odd_v = vld1q_dup_f32(im_x_odd.as_ptr());

        // re * re_odd
        let re_times_re_odd_v = vmulq_f32(re_v, re_odd_v);
        // im * im_odd
        let im_times_im_odd_v = vmulq_f32(im_v, im_odd_v);

        // re * re_odd - im * im_odd
        let re_re_odd_minus_im_im_odd = vsubq_f32(re_times_re_odd_v, im_times_im_odd_v);

        // even_fft[2 * k] + re * re_odd - im * im_odd;
        let result1 = norm(vaddq_f32(
            vld1q_f32(re_x_even.as_ptr()),
            re_re_odd_minus_im_im_odd,
        ));

        out[2 * k] = result1[0];
        out[2 * (k + 1)] = result1[1];
        out[2 * (k + 2)] = result1[2];
        out[2 * (k + 3)] = result1[3];

        // re * im_odd
        let re_times_im_odd_v = vmulq_f32(re_v, im_odd_v);
        // im * re_odd
        let im_times_re_odd_v = vmulq_f32(im_v, re_odd_v);

        // re * im_odd + im * re_odd;
        let re_im_odd_plus_im_re_odd = vaddq_f32(re_times_im_odd_v, im_times_re_odd_v);

        let result2 = norm(vaddq_f32(
            vld1q_f32(im_x_even.as_ptr()),
            re_im_odd_plus_im_re_odd,
        ));

        out[2 * k + 1] = result2[0];
        out[2 * (k + 1) + 1] = result2[1];
        out[2 * (k + 2) + 1] = result2[2];
        out[2 * (k + 3) + 1] = result2[3];

        // out[2 * (k + n_half)] = even_fft[2 * k] - re * re_odd + im * im_odd;
        let result3 = norm(vsubq_f32(
            vld1q_f32(re_x_even.as_ptr()),
            re_re_odd_minus_im_im_odd,
        ));

        out[2 * (k + n_half)] = result3[0];
        out[2 * (k + n_half + 1)] = result3[1];
        out[2 * (k + n_half + 2)] = result3[1];
        out[2 * (k + n_half + 3)] = result3[1];

        // out[2 * (k + n_half) + 1] = even_fft[2 * k + 1] - re * im_odd - im * re_odd;
        let result4 = norm(vsubq_f32(
            vld1q_f32(im_x_even.as_ptr()),
            re_im_odd_plus_im_re_odd,
        ));

        out[2 * (k + n_half) + 1] = result4[0];
        out[2 * (k + n_half + 1) + 1] = result4[1];
        out[2 * (k + n_half + 2) + 1] = result4[1];
        out[2 * (k + n_half + 3) + 1] = result4[1];
    }
}

#[inline(always)]
unsafe fn norm(values: float32x4_t) -> [f32; 4] {
    let mut output = [0.0f32; 4];
    vst1q_f32(output.as_mut_ptr(), values);
    output
}

#[test]
fn fft_test() {
    let m = unsafe {
        fft(&[
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0,
            1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0,
        ])
    };

    // println!("{:?}", m);
}
#[test]
fn tester() {
    unsafe {
        use core::arch::aarch64::*;
        use std::f32::consts::PI;

        let n_t = 40f32;
        let n_quarter = 40;
        let two_pi = PI + PI;
        let pi_vec = vdupq_n_f32(PI);
        let two_pi_vec = vdupq_n_f32(two_pi);
        let n_f32_vec = vdupq_n_f32(n_t);
        let increment = vld1q_f32([0.0, 1.0, 2.0, 3.0].as_ptr());

        for k in (0..n_quarter).step_by(4) {
            let k_base = vdupq_n_f32(k as f32);
            let k_values = vaddq_f32(k_base, increment);
            // pi * [k0..k3]
            let k_pi_values = vmulq_f32(two_pi_vec, k_values);

            // Compute theta_values using vectorized operations
            let s_theta = vdivq_f32(k_pi_values, n_f32_vec);

            // This is a trick such that we can calc cos with sin function. sin(x+pi/2) = cos(x)
            let cos_theta = vaddq_f32(s_theta, vdivq_f32(pi_vec, vdupq_n_f32(2.0)));

            let mut theta_output = [0.0f32; 8];
            vst1q_f32(theta_output.as_mut_ptr(), s_theta);
            vst1q_f32(theta_output.as_mut_ptr().add(4), cos_theta);

            let mut ve = [0.0f32; 8];
            rust_math_neon::sinfv_neon(theta_output.as_mut_ptr(), 8, ve.as_mut_ptr());

            let (_sin, re_) = ve.split_at_mut(4);
            let im_ = vnegq_f32(vld1q_f32(_sin.as_ptr()));
        }
    }
}

// Direct Fourier Transform (DFT) computes the Discrete Fourier Transform directly without any
// optimization. The DFT represents a signal in the frequency domain. It's computationally
// expensive, especially for large inputs, which is why FFT algorithms like Cooley-Tukey are often
// preferred.
fn dft(inp: &[f32]) -> Vec<f32> {
    let zero = 0.0f32;
    let n = inp.len();
    let two_pi = PI + PI;

    let mut out = Vec::new();
    out.reserve(2 * n);
    let n_t = n as f32;

    // Calculate the DFT directly.
    for k in 0..n {
        let k_t = k as f32;
        let mut re = zero;
        let mut im = zero;

        for (j, &inp) in inp.iter().enumerate() {
            let j_t = j as f32;
            let angle = two_pi * k_t * j_t / n_t;
            re += inp * angle.cos();
            im -= inp * angle.sin();
        }

        out.push(re);
        out.push(im);
    }
    out
}

fn main() {}
