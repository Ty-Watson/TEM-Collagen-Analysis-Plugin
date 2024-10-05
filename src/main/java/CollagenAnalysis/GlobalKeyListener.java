package CollagenAnalysis;

import java.awt.KeyboardFocusManager;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;

public class GlobalKeyListener {

    public static boolean esc_pressed = false;
    public static void addGlobalKeyListener() {
        KeyboardFocusManager.getCurrentKeyboardFocusManager().addKeyEventDispatcher(e -> {
            if (e.getID() == KeyEvent.KEY_PRESSED && e.getKeyCode() == KeyEvent.VK_ESCAPE) {
                System.out.println("Escape key pressed globally!");
                esc_pressed = true;
                return true;  // Consume the event
            }
            return false;  // Let other key events proceed
        });
    }
}