
package com.example.multidatasource.config;

import org.springframework.beans.factory.annotation.Qualifier;
import org.springframework.boot.context.properties.ConfigurationProperties;
import org.springframework.boot.jdbc.DataSourceBuilder;
import org.springframework.boot.orm.jpa.EntityManagerFactoryBuilder;
import org.springframework.context.annotation.Bean;
import org.springframework.context.annotation.Configuration;
import org.springframework.context.annotation.Primary;
import org.springframework.data.jpa.repository.config.EnableJpaRepositories;
import org.springframework.orm.jpa.JpaTransactionManager;
import org.springframework.orm.jpa.LocalContainerEntityManagerFactoryBean;
import org.springframework.transaction.PlatformTransactionManager;
import org.springframework.transaction.annotation.EnableTransactionManagement;

import javax.persistence.EntityManagerFactory;
import javax.sql.DataSource;
import java.util.HashMap;
import java.util.Map;

@Configuration
@EnableTransactionManagement
public class DataSourceConfig {

    @Primary
    @Bean(name = "dataSourceA")
    @ConfigurationProperties(prefix = "spring.datasource.a")
    public DataSource dataSourceA() {
        return DataSourceBuilder.create().build();
    }

    @Bean(name = "dataSourceB")
    @ConfigurationProperties(prefix = "spring.datasource.b")
    public DataSource dataSourceB() {
        return DataSourceBuilder.create().build();
    }

    @Bean(name = "dataSourceC")
    @ConfigurationProperties(prefix = "spring.datasource.c")
    public DataSource dataSourceC() {
        return DataSourceBuilder.create().build();
    }

    @Primary
    @Bean(name = "entityManagerFactoryA")
    public LocalContainerEntityManagerFactoryBean entityManagerFactoryA(
            EntityManagerFactoryBuilder builder,
            @Qualifier("dataSourceA") DataSource dataSource) {
        return builder
                .dataSource(dataSource)
                .packages("com.example.multidatasource.model.a")
                .persistenceUnit("a")
                .properties(jpaPropertiesA())
                .build();
    }

    @Bean(name = "entityManagerFactoryB")
    public LocalContainerEntityManagerFactoryBean entityManagerFactoryB(
            EntityManagerFactoryBuilder builder,
            @Qualifier("dataSourceB") DataSource dataSource) {
        return builder
                .dataSource(dataSource)
                .packages("com.example.multidatasource.model.b")
                .persistenceUnit("b")
                .properties(jpaPropertiesB())
                .build();
    }

    @Bean(name = "entityManagerFactoryC")
    public LocalContainerEntityManagerFactoryBean entityManagerFactoryC(
            EntityManagerFactoryBuilder builder,
            @Qualifier("dataSourceC") DataSource dataSource) {
        return builder
                .dataSource(dataSource)
                .packages("com.example.multidatasource.model.c")
                .persistenceUnit("c")
                .properties(jpaPropertiesC())
                .build();
    }

    private Map<String, Object> jpaPropertiesA() {
        Map<String, Object> props = new HashMap<>();
        props.put("hibernate.hbm2ddl.auto", "update");
        props.put("hibernate.dialect", "org.hibernate.dialect.MySQL8Dialect");
        return props;
    }

    private Map<String, Object> jpaPropertiesB() {
        Map<String, Object> props = new HashMap<>();
        props.put("hibernate.hbm2ddl.auto", "update");
        props.put("hibernate.dialect", "org.hibernate.dialect.PostgreSQLDialect");
        return props;
    }

    private Map<String, Object> jpaPropertiesC() {
        Map<String, Object> props = new HashMap<>();
        props.put("hibernate.hbm2ddl.auto", "none");
        props.put("hibernate.dialect", "org.hibernate.dialect.MySQL8Dialect");
        return props;
    }

    @Primary
    @Bean(name = "transactionManagerA")
    public PlatformTransactionManager transactionManagerA(
            @Qualifier("entityManagerFactoryA") EntityManagerFactory entityManagerFactory) {
        return new JpaTransactionManager(entityManagerFactory);
    }

    @Bean(name = "transactionManagerB")
    public PlatformTransactionManager transactionManagerB(
            @Qualifier("entityManagerFactoryB") EntityManagerFactory entityManagerFactory) {
        return new JpaTransactionManager(entityManagerFactory);
    }

    @Bean(name = "transactionManagerC")
    public PlatformTransactionManager transactionManagerC(
            @Qualifier("entityManagerFactoryC") EntityManagerFactory entityManagerFactory) {
        return new JpaTransactionManager(entityManagerFactory);
    }
}
