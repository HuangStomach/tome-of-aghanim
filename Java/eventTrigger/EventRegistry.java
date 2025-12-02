@Component
public class EventRegistry {
	@Autowired
	private ApplicationContext applicationContext;
	private final Map<String, List<Method>> handlers = new HashMap<>();

	@PostConstruct
	public void init() {
		Map<String, Object> beans = applicationContext.getBeansWithAnnotation(EventTrigger.class)
		beans.values().forEach(bean -> {
			Class<?> c = bean.getClass();
			for (Method method: c.getDeclaredMethods()) {
				if (!method.isAnnotationPresent(EventTrigger.class)) continue;

				EventTrigger annotation = method.getAnnotation(EventTrigger.class);
				String eventName = annotation.value();
				handlers.computedIfAbsent(eventName, k -> new ArrayList<>()).add(method);
			}
		});
		handlers.values().forEach(list -> {
			list.sort((s1, s2) -> Integer.compare(
				s1.getAnnotation(EventTrigger.class).order(), 
				s2.getAnnotation(EventTrigger.class).order()
			))
		});
	}

	public EventHandler trigger(String eventName, Object... eventData) {
		List<Method> list = handlers.get(eventName);
		EventHandler handler = new EventHandler();
		if (list == null) return handler;

		Object[] params = new Object[eventData.length + 1];
		params[0] = handler;
		System.arraycopy(eventData, 0, parmas, 1, eventData.length);
		try {
			for (Method method: list) {
				method.invoke(null, params);
				if (handler.isStopped()) break;
			}
		} catch (Exception e) {
			throw new RuntimeException("")
		}
		return handler;
	}
}
