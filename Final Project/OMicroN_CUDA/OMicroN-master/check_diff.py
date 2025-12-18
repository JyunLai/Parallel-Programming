import sys
import math

def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def compare_files(file_serial, file_parallel, tolerance=1e-5):
    print(f"Comparing:\n  1. {file_serial}\n  2. {file_parallel}")
    print(f"Tolerance: {tolerance}")
    
    try:
        with open(file_serial, 'r') as f1, open(file_parallel, 'r') as f2:
            lines1 = f1.readlines()
            lines2 = f2.readlines()
    except FileNotFoundError as e:
        print(f"Error: File not found - {e}")
        return

    if len(lines1) != len(lines2):
        print(f"FAILED: Line counts do not match! ({len(lines1)} vs {len(lines2)})")
        return

    diff_count = 0
    max_diff = 0.0
    
    # 開始逐行比對
    for i, (l1, l2) in enumerate(zip(lines1, lines2)):
        # 跳過註解或空行
        if l1.strip().startswith('#') or not l1.strip():
            continue
            
        # 將每一行切分成數值陣列
        vals1 = [float(v) for v in l1.split() if is_float(v)]
        vals2 = [float(v) for v in l2.split() if is_float(v)]

        if len(vals1) != len(vals2):
            print(f"Line {i+1} format mismatch.")
            diff_count += 1
            continue

        # 比對這一行的每一個數字
        for j, (v1, v2) in enumerate(zip(vals1, vals2)):
            diff = abs(v1 - v2)
            if diff > max_diff:
                max_diff = diff
            
            if diff > tolerance:
                # 有些欄位是 ID (整數)，必須完全一樣；有些是座標或角度，允許誤差
                # 這裡簡單判斷：如果差異很小但在容許範圍外，視為錯誤
                if diff_count < 5: # 只印出前 5 個錯誤以免洗版
                    print(f"Diff at Line {i+1}, Col {j+1}: Serial={v1}, Parallel={v2}, Diff={diff}")
                diff_count += 1

    print("-" * 30)
    if diff_count == 0:
        print(f"SUCCESS: All values match within tolerance {tolerance}.")
        print(f"Max observed difference: {max_diff}")
    else:
        print(f"FAILED: Found {diff_count} discrepancies.")
        print(f"Max observed difference: {max_diff}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python3 check_diff.py <serial_output_file> <parallel_output_file>")
    else:
        compare_files(sys.argv[1], sys.argv[2])