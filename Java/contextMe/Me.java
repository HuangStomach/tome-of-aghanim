import org.springframework.stereotype.Component;
import org.springframework.web.context.annotation.RequestScope;
import lombok.Data;

@Data
@Component
@RequestScope
public class Me {
	private String userId;

	public void destory() {
		// do sth
	}
}