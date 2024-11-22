import java.math.BigInteger;
import java.security.MessageDigest;

public class Hashing {
	public static String GetHash(String s) {
		String result = "";
		
		try {
			MessageDigest md = MessageDigest.getInstance("MD5");
			md.update(s.getBytes());
			result = new BigInteger(1,md.digest()).toString(16);
		} catch (Exception e) {
			e.printStackTrace();
		}
		
		return result;
	}
}
