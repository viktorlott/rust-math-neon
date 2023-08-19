extern crate bindgen;
extern crate cc;

use std::env;
use std::path::PathBuf;

fn main() {
    // // Check if we're compiling for aarch64
    // if let Ok(target) = env::var("TARGET") {
    //     if target != "aarch64-unknown-linux-gnu" {
    //         // Adjust this if your target string is different.
    //         println!("This build script is set up for aarch64 only.");
    //         return;
    //     }
    // } else {
    //     println!("TARGET environment variable not set.");
    //     return;
    // }
    // Compile the C library
    cc::Build::new()
        .files([
            // "math-neon/src/math_acosf.c",
            // "math-neon/src/math_asinf.c",
            // "math-neon/src/math_atan2f.c",
            // "math-neon/src/math_atanf.c",
            // "math-neon/src/math_ceilf.c",
            // "math-neon/src/math_cosf.c",
            // "math-neon/src/math_coshf.c",
            // "math-neon/src/math_expf.c",
            // "math-neon/src/math_fabsf.c",
            // "math-neon/src/math_floorf.c",
            // "math-neon/src/math_fmodf.c",
            // "math-neon/src/math_invsqrtf.c",
            // "math-neon/src/math_ldexpf.c",
            // "math-neon/src/math_log10f.c",
            // "math-neon/src/math_logf.c",
            // "math-neon/src/math_mat2.c",
            // "math-neon/src/math_mat3.c",
            // "math-neon/src/math_mat4.c",
            // "math-neon/src/math_modf.c",
            // "math-neon/src/math_powf.c",
            // "math-neon/src/math_runfast.c",
            // "math-neon/src/math_sincosf.c",
            "math-neon/src/math_sinf.c",
            "math-neon/src/math_sinfv.c",
            // "math-neon/src/math_sinhf.c",
            // "math-neon/src/math_sqrtf.c",
            // "math-neon/src/math_sqrtfv.c",
            // "math-neon/src/math_tanf.c",
            // "math-neon/src/math_tanhf.c",
            // "math-neon/src/math_vec2.c",
            // "math-neon/src/math_vec3.c",
            // "math-neon/src/math_vec4.c",
        ])
        // Ignore warnings when building (Not smart)
        .flag("-Wunused-command-line-argument")
        .flag("-Wunused-variable")
        .flag("-Wno-macro-redefined")
        .flag("-Wno-unused-parameter")
        .flag("-Wno-return-type")
        // NEON optimization flag
        .flag_if_supported("-march=armv8-a")
        .compile("math-neon");

    let bindings = bindgen::Builder::default()
        .header("math-neon/src/math_neon.h")
        .generate()
        .expect("Unable to generate bindings");

    let out_path = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_path.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
