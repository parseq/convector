package deletionsAnalysis;

import java.util.HashMap;

/**
 * Created by german on 11.08.14.
 */

public class DefaultDict<K, V> extends HashMap<K, V> {

    Class<V> cls;

    public DefaultDict(Class cls) {
        this.cls = cls;
    }

    @Override
    public V get(Object key) {
        V returnValue = super.get(key);
        if (returnValue == null) {
            try {
                returnValue = cls.newInstance();
            } catch (Exception e) {
                throw new RuntimeException(e);
            }
            this.put((K) key, returnValue);
        }
        return returnValue;
    }
}