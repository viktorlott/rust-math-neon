fn main() {
    let mut y = [0.0, 0.5, 1.0, 2.0];

    let mut x = [0.0f32; 4];
    unsafe { rust_math_neon::sinfv_neon(y.as_mut_ptr(), 4, x.as_mut_ptr()) };

    //[0.0, 0.4794254, 0.84147096, 0.9092969]
    println!("{:?}", x);

    //[0.0, 0.47942555, 0.841471, 0.9092974]
    println!("{:?}", y.iter().map(|x| x.sin()).collect::<Vec<_>>());
}
