# Personal Notes on EMP-OT

## Structrue
- For each pair of OT messages $(m_0, m_1)$, emp-ot has $\text{LSB}(m_n)=b$
- For COT, the sender has $(m_0, \Delta)$ and the receiver has $(m_{b})$. Note that $m_b$ implies receiver's choise bit $b$. 


## Notes
- extract `block` (`__m128i`) as two `int64` values
    ```c++
    uint64_t pos[2];
    pos[0] = _mm_extract_epi64(r2, 0); // _mm_extract_epi64() for extranction
    pos[1] = _mm_extract_epi64(r2, 1);
    ```
- how to call template fucntion of a templated class, see [here](https://stackoverflow.com/questions/7397934/calling-template-function-within-template-class).
    ```c++
    sender->template send_f2k<OTPre<IO>>(ot, io, i); 
    ```
- 